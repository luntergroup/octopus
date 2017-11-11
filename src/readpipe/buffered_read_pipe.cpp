// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "buffered_read_pipe.hpp"

#include <utility>
#include <limits>
#include <algorithm>
#include <iterator>

#include "utils/mappable_algorithms.hpp"

#include <iostream>

namespace octopus {

BufferedReadPipe::BufferedReadPipe(const ReadPipe& source, Config config)
: BufferedReadPipe {source, config, {}}
{}

BufferedReadPipe::BufferedReadPipe(const ReadPipe& source, Config config, std::vector<GenomicRegion> hints)
: source_ {source}
, config_ {config}
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
    hints_.clear();
}

ReadMap BufferedReadPipe::fetch_reads(const GenomicRegion& region) const
{
    setup_buffer(region);
    return copy_overlapped(buffer_, region);
}

void BufferedReadPipe::hint(std::vector<GenomicRegion> hints) const
{
    hints_.clear();
    for (auto& region : hints) {
        hints_[region.contig_name()].insert(std::move(region));
    }
    for (auto& p : hints_) {
        auto covered_regions = extract_covered_regions(p.second);
        p.second.clear();
        p.second.insert(std::make_move_iterator(std::begin(covered_regions)),
                        std::make_move_iterator(std::end(covered_regions)));
        p.second.shrink_to_fit();
    }
}

// private methods

bool BufferedReadPipe::requires_new_fetch(const GenomicRegion& request) const noexcept
{
    return !buffered_region_ || !contains(*buffered_region_, request);
}

void BufferedReadPipe::setup_buffer(const GenomicRegion& request) const
{
    if (requires_new_fetch(request)) {
        const auto max_region = get_max_fetch_region(request);
        buffered_region_ = source_.get().read_manager().find_covered_subregion(max_region, config_.max_buffer_size);
        buffer_ = source_.get().fetch_reads(*buffered_region_);
    }
}

GenomicRegion BufferedReadPipe::get_max_fetch_region(const GenomicRegion& request) const
{
    const auto default_max_region = get_default_max_fetch_region(request);
    if (hints_.empty() || hints_.count(request.contig_name()) == 0 || hints_.at(request.contig_name()).empty()) {
        return default_max_region;
    } else {
        const auto contained_hints = contained_range(hints_.at(request.contig_name()), default_max_region);
        if (empty(contained_hints)) {
            return request;
        } else {
            if (config_.max_hint_gap) {
                const auto hint_spaces = extract_intervening_regions(contained_hints);
                auto big_space_itr = std::find_if(std::cbegin(hint_spaces), std::cend(hint_spaces),
                                                  [this] (const auto& region) {
                                                      return size(region) > *config_.max_hint_gap;
                                                  });
                if (big_space_itr != std::cend(hint_spaces)) {
                    return left_overhang_region(request, *big_space_itr);
                }
            }
            auto hint_region = encompassing_region(contained_hints);
            return closed_region(request, hint_region);
        }
    }
}

namespace {

GenomicRegion fully_expand_rhs(const GenomicRegion& region)
{
    return GenomicRegion {region.contig_name(), region.begin(), std::numeric_limits<GenomicRegion::Position>::max()};
}

} // namespace

GenomicRegion BufferedReadPipe::get_default_max_fetch_region(const GenomicRegion& request) const
{
    if (config_.max_fetch_size) {
        return expand_rhs(request, *config_.max_fetch_size);
    } else {
        return fully_expand_rhs(request);
    }
}

} // namespace octopus
