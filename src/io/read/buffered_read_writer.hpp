// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef buffered_read_writer_hpp
#define buffered_read_writer_hpp

#include <deque>
#include <functional>
#include <cstddef>
#include <utility>
#include <algorithm>
#include <iterator>
#include <cassert>

#include <boost/optional.hpp>

#include "basics/aligned_read.hpp"
#include "utils/memory_footprint.hpp"
#include "utils/append.hpp"
#include "read_writer.hpp"

namespace octopus { namespace io {

template <typename Read>
class BufferedReadWriter
{
public:
    struct Config
    {
        boost::optional<std::size_t> max_buffer_size = boost::none;
        boost::optional<MemoryFootprint> max_buffer_footprint = boost::none;
        enum class WriteOrder { given, region_sorted, bam_sorted } write_order = WriteOrder::bam_sorted;
    };
    
    BufferedReadWriter() = delete;
    
    BufferedReadWriter(ReadWriter& writer);
    BufferedReadWriter(ReadWriter& writer, Config config);
    
    BufferedReadWriter(const BufferedReadWriter&)            = default;
    BufferedReadWriter& operator=(const BufferedReadWriter&) = default;
    BufferedReadWriter(BufferedReadWriter&&)                 = default;
    BufferedReadWriter& operator=(BufferedReadWriter&&)      = default;
    
    ~BufferedReadWriter();
    
    const ReadWriter& writer() const noexcept;
    ReadWriter& writer() noexcept;
    
    Config config() const noexcept;
    
    std::size_t buffer_size() const noexcept;
    MemoryFootprint buffer_footprint() const noexcept;
    
    void write(Read read);
    template <typename Range> void write(Range reads);
    
    void flush();
    void clear() noexcept;

private:
    std::deque<Read> buffer_;
    std::reference_wrapper<ReadWriter> writer_;
    Config config_;
    MemoryFootprint buffer_footprint_;
    
    void sort_buffer();
    bool can_buffer(MemoryFootprint footprint, std::size_t reads = 1) const noexcept;
    bool can_buffer(const Read& read) const noexcept;
    template <typename Range> bool can_buffer_all(const Range& reads) const noexcept;
    void buffer(Read read);
    template <typename Range> void buffer(Range reads);
};

template <typename Read>
BufferedReadWriter<Read>::BufferedReadWriter(ReadWriter& writer) : BufferedReadWriter {writer, Config {}} {};

template <typename Read>
BufferedReadWriter<Read>::BufferedReadWriter(ReadWriter& writer, Config config)
: buffer_ {}
, writer_ {writer}
, config_ {config}
, buffer_footprint_ {0}
{}

template <typename Read>
BufferedReadWriter<Read>::~BufferedReadWriter()
{
    try {
        flush();
    } catch (...) {
        // TODO: log?
        clear();
    }
}

template <typename Read>
const ReadWriter& BufferedReadWriter<Read>::writer() const noexcept
{
    return writer_.get();
}

template <typename Read>
ReadWriter& BufferedReadWriter<Read>::writer() noexcept
{
    return writer_.get();
}

template <typename Read>
typename BufferedReadWriter<Read>::Config BufferedReadWriter<Read>::config() const noexcept
{
    return config_;
}

template <typename Read>
std::size_t BufferedReadWriter<Read>::buffer_size() const noexcept
{
    return buffer_.size();
}

template <typename Read>
MemoryFootprint BufferedReadWriter<Read>::buffer_footprint() const noexcept
{
    return buffer_footprint_;
}

template <typename Read>
void BufferedReadWriter<Read>::write(Read read)
{
    if (can_buffer(read)) {
        buffer(std::move(read));
    } else if (buffer_.empty()) {
        writer_.get().write(read);
    } else {
        flush();
        write(std::move(read));
    }
}

template <typename Read>
void BufferedReadWriter<Read>::flush()
{
    sort_buffer();
    writer_.get() << buffer_;
    buffer_.clear();
    buffer_footprint_ = 0;
}

template <typename Read>
void BufferedReadWriter<Read>::clear() noexcept
{
    buffer_.clear();
    buffer_footprint_ = 0;
}

// private methods

namespace  {

template <typename Container>
void sort(Container& reads)
{
    std::sort(std::begin(reads), std::end(reads));
}

struct BAMLess
{
    bool operator()(const AlignedRead& lhs, const AlignedRead& rhs) const noexcept
    {
        return clipped_mapped_region(lhs) < clipped_mapped_region(rhs);
    }
};

template <typename Container>
void bam_sort(Container& reads)
{
    std::sort(std::begin(reads), std::end(reads), BAMLess {});
}

} // namespace

template <typename Read>
void BufferedReadWriter<Read>::sort_buffer()
{
    switch (config_.write_order) {
        case Config::WriteOrder::region_sorted: sort(buffer_); break;
        case Config::WriteOrder::bam_sorted: bam_sort(buffer_); break;
        default: break;
    }
}

template <typename Read>
bool BufferedReadWriter<Read>::can_buffer(const MemoryFootprint footprint, const std::size_t reads) const noexcept
{
    return (!config_.max_buffer_size || (buffer_size() + reads) <= config_.max_buffer_size)
           && (!config_.max_buffer_footprint || (buffer_footprint() + footprint) <= config_.max_buffer_footprint);
}

template <typename Read>
bool BufferedReadWriter<Read>::can_buffer(const Read& read) const noexcept
{
    return (buffer_.empty() || is_same_contig(buffer_.front(), read)) && can_buffer(footprint(read));
}

template <typename Read>
void BufferedReadWriter<Read>::buffer(Read read)
{
    buffer_.push_back(std::move(read));
    buffer_footprint_ += footprint(buffer_.back());
}

// non-member methods

template <typename Read>
BufferedReadWriter<Read>& operator<<(BufferedReadWriter<Read>& dst, const Read& read)
{
    dst.write(read);
    return dst;
}

template <typename Read>
template <typename Range>
void BufferedReadWriter<Read>::write(Range reads)
{
    if (reads.empty()) return;
    if (can_buffer_all(reads)) {
        buffer(std::move(reads));
    } else if (!buffer_.empty() && !is_same_contig(buffer_.front(), reads.front())) {
        flush();
        write(std::move(reads));
    } else {
        utils::append(std::move(reads), buffer_);
        sort_buffer();
        auto buffer_flush_end_itr = std::cbegin(buffer_);
        for (auto buffer_itr = std::cbegin(buffer_); buffer_itr != std::cend(buffer_);) {
            buffer_footprint_ = 0;
            buffer_itr = std::find_if(buffer_itr, std::cend(buffer_), [this] (const auto& read) {
                if (can_buffer(read)) {
                    buffer_footprint_ += footprint(read);
                    return false;
                } else {
                    return true;
                }
            });
            if (buffer_itr != std::cend(buffer_)) {
                buffer_flush_end_itr = buffer_itr;
            }
        }
        std::for_each(std::cbegin(buffer_), buffer_flush_end_itr, [this] (const auto& read) { writer_.get().write(read); });
        buffer_.erase(std::cbegin(buffer_), buffer_flush_end_itr);
        buffer_footprint_ = footprint(buffer_);
    }
}

template <typename Read>
template <typename Range>
bool BufferedReadWriter<Read>::can_buffer_all(const Range& reads) const noexcept
{
    return reads.empty() || ((buffer_.empty() || is_same_contig(buffer_.front(), reads.front())) && can_buffer(footprint(reads), reads.size()));
}

template <typename Read>
template <typename Range>
void BufferedReadWriter<Read>::buffer(Range reads)
{
    buffer_footprint_ += footprint(reads);
    utils::append(std::move(reads), buffer_);
}

template <typename Read, typename Range>
BufferedReadWriter<Read>& operator<<(BufferedReadWriter<Read>& dst, const Range& reads)
{
    dst.write(reads);
    return dst;
}

} // namespace io
} // namespace octopus

#endif
