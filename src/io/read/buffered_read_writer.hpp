// Copyright (c) 2017 Daniel Cooke
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
    
    void write(AlignedRead read);
    template <typename Range> void write(Range reads);
    
    void flush();
    void clear() noexcept;

private:
    std::deque<AlignedRead> buffer_;
    std::reference_wrapper<ReadWriter> writer_;
    Config config_;
    MemoryFootprint buffer_footprint_;
    
    void sort_buffer();
    bool can_buffer(MemoryFootprint footprint, std::size_t reads = 1) const noexcept;
    bool can_buffer(const AlignedRead& read) const noexcept;
    template <typename Range> bool can_buffer(const Range& reads) const noexcept;
    void buffer(AlignedRead read);
    template <typename Range> void buffer(Range reads);
};

template <typename Range>
void BufferedReadWriter::write(Range reads)
{
    if (reads.empty()) return;
    if (can_buffer(reads)) {
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
            buffer_itr = std::find_if(buffer_itr, std::cend(buffer_), [this] (const AlignedRead& read) {
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
        std::for_each(std::cbegin(buffer_), buffer_flush_end_itr, [this] (const AlignedRead& read) { writer_.get().write(read); });
        buffer_.erase(std::cbegin(buffer_), buffer_flush_end_itr);
        buffer_footprint_ = footprint(buffer_);
    }
}

template <typename Range>
bool BufferedReadWriter::can_buffer(const Range& reads) const noexcept
{
    assert(!reads.empty());
    return (buffer_.empty() || is_same_contig(buffer_.front(), reads.front())) && can_buffer(footprint(reads), reads.size());
}

template <typename Range>
void BufferedReadWriter::buffer(Range reads)
{
    buffer_footprint_ += footprint(reads);
    utils::append(std::move(reads), buffer_);
}

BufferedReadWriter& operator<<(BufferedReadWriter& dst, const AlignedRead& read);

template <typename Range>
BufferedReadWriter& operator<<(BufferedReadWriter& dst, const Range& reads)
{
    dst.write(reads);
    return dst;
}

} // namespace io
} // namespace octopus

#endif
