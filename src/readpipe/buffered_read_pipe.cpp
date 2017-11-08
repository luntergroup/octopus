// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "buffered_read_pipe.hpp"

#include "utils/mappable_algorithms.hpp"

namespace octopus {

BufferedReadPipe::BufferedReadPipe(const ReadPipe& source, std::size_t max_buffer_size)
: source_ {source}
, max_buffer_size_ {max_buffer_size}
, buffer_ {}
, buffered_region_ {}
, hints_ {}
{}

BufferedReadPipe::BufferedReadPipe(const ReadPipe& source, std::size_t max_buffer_size,
                                   std::vector<GenomicRegion> hints)
: source_ {source}
, max_buffer_size_ {max_buffer_size}
, buffer_ {}
, buffered_region_ {}
, hints_ {}
{
    hint(std::move(hints));
}

const ReadPipe& BufferedReadPipe::source() const noexcept
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

void BufferedReadPipe::hint(std::vector<GenomicRegion> hints) const
{
    for (auto& region : hints) {
        hints_[region.contig_name()].insert(std::move(region));
    }
}

// private methods

bool BufferedReadPipe::requires_new_fetch(const GenomicRegion& region) const noexcept
{
    return !buffered_region_ || !contains(*buffered_region_, region);
}

void BufferedReadPipe::setup_buffer(const GenomicRegion& region) const
{
    if (requires_new_fetch(region)) {
        const auto max_region = calculate_max_buffer_region(region);
        buffered_region_ = source_.get().read_manager().find_covered_subregion(max_region, max_buffer_size_);
        buffer_ = source_.get().fetch_reads(*buffered_region_);
    }
}

GenomicRegion BufferedReadPipe::calculate_max_buffer_region(const GenomicRegion& request) const
{
    const auto default_max_region = expand_rhs(request, 2 * max_buffer_size_);
    if (hints_.empty() || hints_.at(request.contig_name()).empty()) {
        return default_max_region;
    } else {
        const auto contained_hints = contained_range(hints_.at(request.contig_name()), default_max_region);
        if (empty(contained_hints)) {
            return request;
        } else {
            auto hint_region = encompassing_region(contained_hints);
            return closed_region(request, hint_region);
        }
    }
}

} // namespace octopus
