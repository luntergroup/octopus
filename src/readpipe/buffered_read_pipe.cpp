// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "buffered_read_pipe.hpp"

#include <utility>
#include <limits>
#include <algorithm>
#include <iterator>

#include "utils/mappable_algorithms.hpp"
#include "utils/read_stats.hpp"

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
, debug_log_ {}
{
    hint(std::move(hints));
    if (DEBUG_MODE) debug_log_ = logging::DebugLogger {};
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
    if (config_.max_buffer_size == 0) return source_.get().fetch_reads(region);
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

bool BufferedReadPipe::is_cached(const GenomicRegion& region) const noexcept
{
    return buffered_region_ && contains(*buffered_region_, region);
}

// private methods

void BufferedReadPipe::setup_buffer(const GenomicRegion& request) const
{
    if (!is_cached(request)) {
        if (debug_log_) stream(*debug_log_) << "Request " << request << " is not cached";
        auto max_region = get_max_fetch_region(request);
        if (debug_log_) stream(*debug_log_) << "Max fetch region for request " << request << " is " << max_region;
        bool unchecked_fetch {false};
        if (can_make_unchecked_fetch()) {
            buffered_region_ = std::move(max_region);
            unchecked_fetch = true;
        } else {
            buffered_region_ = source_.get().read_manager().find_covered_subregion(max_region, config_.max_buffer_size);
        }
        if (debug_log_) stream(*debug_log_) << "Buffer region for request " << request << " is " << *buffered_region_;
        buffer_ = source_.get().fetch_reads(expand(*buffered_region_, config_.fetch_expansion));
        if (unchecked_fetch) {
            const auto fetch_size = count_reads(buffer_);
            if (fetch_size > config_.max_buffer_size) {
                if (default_unchecked_fetch_overflowed_) {
                    adjusted_unchecked_fetch_overflowed_ = true;
                } else {
                    default_unchecked_fetch_overflowed_ = true;
                }
                // Clear buffer of reads to rhs of request
                for (auto& p : buffer_) {
                    const auto last_overlapped = find_first_after(p.second, request);
                    p.second.erase(last_overlapped, std::cend(p.second));
                }
                buffered_region_ = request;
            }
        } else {
            if (min_checked_fetch_size_) {
                min_checked_fetch_size_ = std::min(size(*buffered_region_), *min_checked_fetch_size_);
            } else {
                min_checked_fetch_size_ = size(*buffered_region_);
            }
        }
    } else if (debug_log_) {
        stream(*debug_log_) << "Request " << request << " is already cached";
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
    } else if (min_checked_fetch_size_) {
        if (size(request) < *min_checked_fetch_size_) {
            const auto recommended_expansion = *min_checked_fetch_size_ - size(request);
            return expand_rhs(request, recommended_expansion);
        } else {
            return request;
        }
    } else {
        return fully_expand_rhs(request);
    }
}

bool BufferedReadPipe::can_make_unchecked_fetch() const noexcept
{
    return config_.allow_unchecked_fetches
           && !(default_unchecked_fetch_overflowed_ && !min_checked_fetch_size_)
           && !adjusted_unchecked_fetch_overflowed_
           && (config_.max_fetch_size || min_checked_fetch_size_);
}

} // namespace octopus
