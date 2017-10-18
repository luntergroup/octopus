// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "buffered_read_pipe.hpp"

namespace octopus {

BufferedReadPipe::BufferedReadPipe(ReadPipe& source, std::size_t max_buffer_size)
: source_ {source}
, max_buffer_size_ {max_buffer_size}
, buffer_ {}
, buffered_region_ {}
{}

ReadPipe& BufferedReadPipe::source() noexcept
{
    return source_.get();
}

void BufferedReadPipe::clear() noexcept
{
    buffer_.clear();
    buffered_region_ = boost::none;
}

ReadMap BufferedReadPipe::fetch_reads(const GenomicRegion& region) const
{
    setup_buffer(region);
    return copy_overlapped(buffer_, region);
}

// private methods

bool BufferedReadPipe::requires_new_fetch(const GenomicRegion& region) const noexcept
{
    return !buffered_region_ || !contains(*buffered_region_, region);
}

void BufferedReadPipe::setup_buffer(const GenomicRegion& region) const
{
    if (requires_new_fetch(region)) {
        const auto max_region = expand_rhs(region, 2 * max_buffer_size_);
        buffered_region_ = source_.get().read_manager().find_covered_subregion(max_region, max_buffer_size_);
        buffer_ = source_.get().fetch_reads(*buffered_region_);
    }
}

} // namespace octopus
