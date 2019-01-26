// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "buffered_read_writer.hpp"

namespace octopus { namespace io {

// public methods

BufferedReadWriter::BufferedReadWriter(ReadWriter& writer) : BufferedReadWriter {writer, Config {}} {};

BufferedReadWriter::BufferedReadWriter(ReadWriter& writer, Config config)
: buffer_ {}
, writer_ {writer}
, config_ {config}
, buffer_footprint_ {0}
{}

BufferedReadWriter::~BufferedReadWriter()
{
    try {
        flush();
    } catch (...) {
        // TODO: log?
        clear();
    }
}

const ReadWriter& BufferedReadWriter::writer() const noexcept
{
    return writer_.get();
}

ReadWriter& BufferedReadWriter::writer() noexcept
{
    return writer_.get();
}

BufferedReadWriter::Config BufferedReadWriter::config() const noexcept
{
    return config_;
}

std::size_t BufferedReadWriter::buffer_size() const noexcept
{
    return buffer_.size();
}

MemoryFootprint BufferedReadWriter::buffer_footprint() const noexcept
{
    return buffer_footprint_;
}

void BufferedReadWriter::write(AlignedRead read)
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

void BufferedReadWriter::flush()
{
    sort_buffer();
    writer_.get() << buffer_;
    buffer_.clear();
    buffer_footprint_ = 0;
}

void BufferedReadWriter::clear() noexcept
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

void BufferedReadWriter::sort_buffer()
{
    switch (config_.write_order) {
        case Config::WriteOrder::region_sorted: sort(buffer_); break;
        case Config::WriteOrder::bam_sorted: bam_sort(buffer_); break;
        default: break;
    }
}

bool BufferedReadWriter::can_buffer(const MemoryFootprint footprint, const std::size_t reads) const noexcept
{
    return (!config_.max_buffer_size || (buffer_size() + reads) <= config_.max_buffer_size)
           && (!config_.max_buffer_footprint || (buffer_footprint() + footprint) <= config_.max_buffer_footprint);
}

bool BufferedReadWriter::can_buffer(const AlignedRead& read) const noexcept
{
    return (buffer_.empty() || is_same_contig(buffer_.front(), read)) && can_buffer(footprint(read));
}

void BufferedReadWriter::buffer(AlignedRead read)
{
    buffer_.push_back(std::move(read));
    buffer_footprint_ += footprint(buffer_.back());
}

// non-member methods

BufferedReadWriter& operator<<(BufferedReadWriter& dst, const AlignedRead& read)
{
    dst.write(read);
    return dst;
}

} // namespace io
} // namespace octopus
