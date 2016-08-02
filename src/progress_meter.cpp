//
//  progress_meter.cpp
//  Octopus
//
//  Created by Daniel Cooke on 10/03/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#include "progress_meter.hpp"

#include <sstream>
#include <algorithm>
#include <functional>
#include <cmath>
#include <utility>
#include <initializer_list>
#include <cassert>

#include "mappable_algorithms.hpp"
#include "mappable_map.hpp"
#include "timing.hpp"
#include "string_utils.hpp"
#include "maths.hpp"

namespace octopus {

using utils::TimeInterval;

namespace
{
    template <typename T>
    unsigned num_digits(T x)
    {
        return static_cast<unsigned>(std::to_string(x).size());
    }
} // namespace

template <typename T>
unsigned max_str_length(const T& p)
{
    auto result = static_cast<unsigned>(p.first.size());
    if (p.second.empty()) return result;
    return result + num_digits(p.second.rightmost().end());
}

unsigned max_position_str_length(const InputRegionMap& input_regions)
{
    assert(!input_regions.empty());
    
    unsigned result {0};
    
    for (const auto& p : input_regions) {
        const auto curr = max_str_length(p);
        if (result < curr) result = curr;
    }
    
    return result;
}

unsigned max_position_str_length(const GenomicRegion& region)
{
    return static_cast<unsigned>(region.contig_name().size()) + num_digits(region.end());
}

template <typename T>
unsigned calculate_position_tab_length(const T& region)
{
    return std::max(18u, max_position_str_length(region));
}

ProgressMeter::ProgressMeter(InputRegionMap regions)
:
regions_ {std::move(regions)},
completed_regions_ {},
num_bp_to_search_ {sum_region_sizes(regions_)},
num_bp_completed_ {0},
percent_unitl_log_ {percent_block_size_},
percent_at_last_log_ {0},
start_ {std::chrono::system_clock::now()},
last_log_ {start_},
done_ {false},
position_tab_length_ {},
block_compute_times_ {},
log_ {}
{
    if (!regions_.empty()) {
        position_tab_length_ = calculate_position_tab_length(regions_);
    }
}

ProgressMeter::ProgressMeter(GenomicRegion region)
:
ProgressMeter {
    InputRegionMap {std::make_pair(region.contig_name(), InputRegionMap::mapped_type {region})}
}
{}

ProgressMeter::ProgressMeter(ProgressMeter&& other)
{
    using std::move;
    std::lock_guard<std::mutex> lock {other.mutex_};
    regions_ = move(other.regions_);
    completed_regions_ = move(other.completed_regions_);
    num_bp_to_search_ = move(other.num_bp_to_search_);
    num_bp_completed_ = move(other.num_bp_completed_);
    percent_block_size_ = move(other.percent_block_size_);
    percent_unitl_log_ = move(other.percent_unitl_log_);
    percent_at_last_log_ = move(other.percent_at_last_log_);
    start_ = move(other.start_);
    last_log_ = move(other.last_log_);
    done_ = move(other.done_);
    position_tab_length_ = move(other.position_tab_length_);
    block_compute_times_ = move(other.block_compute_times_);
    log_ = move(other.log_);
}

double percent_completed(const std::size_t num_bp_completed, const std::size_t num_bp_to_search)
{
    return 100 * static_cast<double>(num_bp_completed) / num_bp_to_search;
}

auto percent_completed_str(const std::size_t num_bp_completed, const std::size_t num_bp_to_search)
{
    return utils::to_string(percent_completed(num_bp_completed, num_bp_to_search), 1) + '%';
}

std::string to_string(const TimeInterval& duration)
{
    std::ostringstream ss {};
    ss << duration;
    return ss.str();
}

template <typename ForwardIt>
auto mean_duration(ForwardIt first, ForwardIt last)
{
    return maths::mean(first, last, [] (const auto& d) { return d.count(); });
}

template <typename Container>
auto mean_duration(const Container& durations)
{
    return mean_duration(std::cbegin(durations), std::cend(durations));
}

template <typename ForwardIt>
auto stdev_duration(ForwardIt first, ForwardIt last)
{
    return maths::stdev(first, last, [] (const auto& d) { return d.count(); });
}

template <typename Container>
auto stdev_duration(const Container& durations)
{
    return stdev_duration(std::cbegin(durations), std::cend(durations));
}

template <typename Container>
void remove_outliers(Container& durations)
{
    if (durations.size() < 2) {
        return;
    }
    
    if (std::adjacent_find(std::cbegin(durations), std::cend(durations),
        std::not_equal_to<void> {}) == std::cend(durations)) {
        return;
    }
    
    auto it = std::min_element(std::begin(durations), std::end(durations));
    
    if (it == std::begin(durations)) {
        it = std::remove(it, std::end(durations), *it);
    } else {
        it = std::end(durations);
    }
    
    const auto mean  = mean_duration(std::begin(durations), it);
    const auto stdev = stdev_duration(std::begin(durations), it);
    
    const auto min = std::max(0.0, mean - (2 * stdev));
    const auto max = mean + (2 * stdev);
    
    it = std::remove_if(it, std::end(durations),
                        [=] (const auto& duration) {
                            return duration.count() < min || duration.count() > max;
                        });
    
    durations.erase(it, std::end(durations));
}

template <typename Container>
auto estimate_ttc(const std::chrono::system_clock::time_point now,
                  const Container& durations,
                  const std::size_t num_remaining_blocks)
{
    const auto mean_block_duration = mean_duration(durations);
    
    const auto estimated_remaining_ticks = static_cast<std::size_t>(num_remaining_blocks * mean_block_duration);
    
    const std::chrono::milliseconds estimated_remaining_duration {estimated_remaining_ticks};
    
    return TimeInterval {now, now + estimated_remaining_duration};
}

ProgressMeter::~ProgressMeter()
{
    if (!done_ && !regions_.empty() && num_bp_completed_ > 0) {
        const TimeInterval duration {start_, std::chrono::system_clock::now()};
        
        const auto time_taken = to_string(duration);
        
        stream(log_) << std::string(position_tab_length_ - 4, ' ')
                     << "-"
                     << completed_pad("100%")
                     << "100%"
                     << time_taken_pad(time_taken)
                     << time_taken
                     << ttc_pad("-")
                     << "-";
    }
}

void ProgressMeter::set_percent_block_size(const double size)
{
    block_compute_times_.clear(); // TODO: can we use old block times to estimate new block times?
    percent_block_size_ = size;
    percent_unitl_log_  = std::min(size, percent_unitl_log_);
}

void ProgressMeter::start()
{
    if (!regions_.empty()) {
        completed_regions_.reserve(regions_.size());
        write_header();
    }
    start_ = std::chrono::system_clock::now();
    last_log_ = start_;
}

void ProgressMeter::resume()
{
    // TODO
}

void ProgressMeter::pause()
{
    // TODO
}

void ProgressMeter::stop()
{
    if (!done_ && !regions_.empty()) {
        const TimeInterval duration {start_, std::chrono::system_clock::now()};
        
        const auto time_taken = to_string(duration);
        
        stream(log_) << std::string(position_tab_length_ - 4, ' ')
        << "-"
        << completed_pad("100%")
        << "100%"
        << time_taken_pad(time_taken)
        << time_taken
        << ttc_pad("-")
        << "-";
    }
    
    done_ = true;
}

void ProgressMeter::log_completed(const GenomicRegion& region)
{
    std::lock_guard<std::mutex> lock {mutex_};
    
    const auto new_bp_processed = merge(region);
    
    const auto new_percent_done = percent_completed(new_bp_processed, num_bp_to_search_);
    
    num_bp_completed_ += new_bp_processed;
    
    percent_unitl_log_ -= new_percent_done;
    
    if (percent_unitl_log_ <= 0) {
        output_log(region);
    }
}

// private methods

ProgressMeter::RegionSizeType ProgressMeter::merge(const GenomicRegion& region)
{
    const auto& contig_region = region.contig_region();
    
    if (completed_regions_.count(region.contig_name()) == 0) {
        completed_regions_.emplace(std::piecewise_construct,
                                   std::forward_as_tuple(region.contig_name()),
                                   std::forward_as_tuple(std::initializer_list<ContigRegion>({contig_region})));
        return size(contig_region);
    }
    
    auto& completed_regions = completed_regions_.at(region.contig_name());
    
    if (completed_regions.has_overlapped(contig_region)) {
        auto overlapped = completed_regions.overlap_range(contig_region);
        
        if (size(overlapped) == 1) {
            RegionSizeType result {0};
            
            if (contains(overlapped.front(), contig_region)) return result;
            
            auto new_region = contig_region;
            
            if (begins_before(new_region, overlapped.front())) {
                result += left_overhang_size(new_region, overlapped.front());
            } else if (begins_before(overlapped.front(), new_region)) {
                new_region = encompassing_region(overlapped.front(), new_region);
            }
            
            if (ends_before(overlapped.front(), new_region)) {
                result += right_overhang_size(new_region, overlapped.front());
            } else if (ends_before(new_region, overlapped.front())) {
                new_region = encompassing_region(overlapped.front(), new_region);
            }
            
            completed_regions.erase(overlapped.front());
            completed_regions.insert(std::move(new_region));
            
            return result;
        }
        
        auto new_region = contig_region;
        auto result     = size(new_region);
        
        result -= overlap_size(overlapped.front(), contig_region);
        
        if (begins_before(overlapped.front(), new_region)) {
            new_region = encompassing_region(overlapped.front(), new_region);
        }
        
        overlapped.advance_begin(1);
        
        result -= overlap_size(overlapped.back(), new_region);
        
        if (ends_before(new_region, overlapped.back())) {
            new_region = encompassing_region(new_region, overlapped.back());
        }
        
        overlapped.advance_end(-1);
        
        for (const auto& r : overlapped) result -= size(r);
        
        completed_regions.erase_overlapped(contig_region);
        completed_regions.insert(std::move(new_region));
        
        return result;
    }
    
    completed_regions.insert(contig_region);
    
    return size(contig_region);
}

void ProgressMeter::write_header()
{
    const std::string pos_tab_bar(position_tab_length_, '-');
    
    const auto num_position_tab_whitespaces = position_tab_length_ - 8;
    
    const std::string pos_tab_lhs_pad((position_tab_length_ - 8) / 2, ' ');
    
    auto pos_tab_rhs_pad = pos_tab_lhs_pad;
    
    if (num_position_tab_whitespaces % 2) {
        pos_tab_rhs_pad += ' ';
    }
    
    stream(log_) << pos_tab_bar << "------------------------------------------------------";
    stream(log_) << pos_tab_lhs_pad << "current " << pos_tab_rhs_pad
    << "|                   |     time      |     estimated   ";
    stream(log_) << pos_tab_lhs_pad << "position" << pos_tab_rhs_pad
    << "|     completed     |     taken     |     ttc         ";
    stream(log_) << pos_tab_bar << "------------------------------------------------------";
}

void ProgressMeter::output_log(const GenomicRegion& region)
{
    const auto percent_done = percent_completed(num_bp_completed_, num_bp_to_search_);
    
    const auto now = std::chrono::system_clock::now();
    
    const TimeInterval duration {start_, now};
    
    const auto time_taken = to_string(duration);
    
    if (percent_done >= 100) return;
    
    const auto percent_since_last_log = percent_done - percent_at_last_log_;
    
    const auto num_blocks_completed = static_cast<std::size_t>(std::floor(percent_since_last_log / percent_block_size_));
    
    std::string ttc {"-"};
    
    if (num_blocks_completed > 0) {
        const auto duration_since_last_log = std::chrono::duration_cast<DurationUnits>(now - last_log_);
        
        const DurationUnits duration_per_block {duration_since_last_log.count() / num_blocks_completed};
        
        std::fill_n(std::back_inserter(block_compute_times_), num_blocks_completed,
                    duration_per_block);
        
        const auto num_remaining_blocks = static_cast<std::size_t>((100.0 - percent_done) / percent_block_size_);
        
        remove_outliers(block_compute_times_);
        
        ttc = to_string(estimate_ttc(now, block_compute_times_, num_remaining_blocks));
        
        assert(!ttc.empty());
        
        if (ttc.front() == '0') {
            ttc = "-";
        }
    }
    
    const auto percent_completed = percent_completed_str(num_bp_completed_, num_bp_to_search_);
    
    stream(log_) << curr_position_pad(region)
                 << region.contig_name() << ':' << region.end()
                 << completed_pad(percent_completed)
                 << percent_completed
                 << time_taken_pad(time_taken)
                 << time_taken
                 << ttc_pad(ttc)
                 << ttc;
    
    last_log_            = now;
    percent_unitl_log_   = percent_block_size_;
    percent_at_last_log_ = percent_done;
}

std::string ProgressMeter::curr_position_pad(const GenomicRegion& completed_region) const
{
    assert(position_tab_length_ > 3);
    
    const auto num_contig_name_letters = completed_region.contig_name().size();
    const auto num_region_end_digits = num_digits(completed_region.end());
    
    const auto l = num_contig_name_letters + num_region_end_digits + 1; // +1 for ':'
    
    if (l < position_tab_length_ - 3) {
        return std::string(position_tab_length_ - l - 3, ' ');
    }
    
    return "";
}

std::string ProgressMeter::completed_pad(const std::string& percent_completed) const
{
    return std::string(std::size_t {17} - percent_completed.size(), ' ');
}

std::string ProgressMeter::time_taken_pad(const std::string& time_taken) const
{
    return std::string(16 - time_taken.size(), ' ');
}

std::string ProgressMeter::ttc_pad(const std::string& ttc) const
{
    return std::string(16 - ttc.size(), ' ');
}

} // namespace octopus
