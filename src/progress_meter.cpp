//
//  progress_meter.cpp
//  Octopus
//
//  Created by Daniel Cooke on 10/03/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#include "progress_meter.hpp"

#include <string>
#include <sstream>
#include <algorithm>
#include <functional>
#include <cmath>
#include <cassert>

#include "mappable_algorithms.hpp"
#include "mappable_map.hpp"
#include "timing.hpp"
#include "string_utils.hpp"
#include "maths.hpp"

namespace Octopus
{
    template <typename L>
    void write_header(L& log)
    {
        log << "------------------------------------------------------------------------";
        log << "     current      |                   |     time      |     estimated   ";
        log << "     position     |     completed     |     taken     |     ttc         ";
        log << "------------------------------------------------------------------------";
    }
    
    ProgressMeter::ProgressMeter(const SearchRegions& input_regions)
    :
    num_bp_completed_ {0},
    start_ {std::chrono::system_clock::now()},
    last_log_ {std::chrono::system_clock::now()}
    {
        write_header(log_);
    }
    
    ProgressMeter::ProgressMeter(const GenomicRegion& input_region)
    :
    num_bp_completed_ {0},
    percent_unitl_log_ {percent_block_size_},
    percent_at_last_log_ {0},
    start_ {std::chrono::system_clock::now()},
    last_log_ {std::chrono::system_clock::now()},
    done_ {false}
    {
        regions_[input_region.get_contig_name()].emplace(input_region);
        
        completed_region_ = head_region(input_region);
        
        num_bp_to_search_ = sum_region_sizes(regions_);
        
        write_header(log_);
    }
    
    template <typename T>
    auto num_digits(T x)
    {
        return std::to_string(x).size();
    }
    
    std::string curr_position_pad(const GenomicRegion& completed_region)
    {
        const auto num_contig_name_letters = completed_region.get_contig_name().size();
        const auto num_region_begin_digits = num_digits(completed_region.get_begin());
        if (num_contig_name_letters + num_region_begin_digits <= 15) {
            return std::string(15 - (num_contig_name_letters + num_region_begin_digits), ' ');
        } else {
            return "";
        }
    }
    
    double percent_completed(const std::size_t num_bp_completed,
                             const std::size_t num_bp_to_search)
    {
        return 100 * static_cast<double>(num_bp_completed) / num_bp_to_search;
    }
    
    auto percent_completed_str(const std::size_t num_bp_completed,
                               const std::size_t num_bp_to_search)
    {
        return Octopus::to_string(percent_completed(num_bp_completed, num_bp_to_search), 1) + '%';
    }
    
    std::string completed_pad(const std::string& percent_completed)
    {
        return std::string(std::size_t {17} - percent_completed.size(), ' ');
    }
    
    std::string to_string(const TimeInterval& duration)
    {
        std::ostringstream ss {};
        ss << duration;
        return ss.str();
    }
    
    std::string time_taken_pad(const std::string& time_taken)
    {
        return std::string(16 - time_taken.size(), ' ');
    }
    
    template <typename ForwardIt>
    auto mean_duration(ForwardIt first, ForwardIt last)
    {
        return Maths::mean(first, last, [] (const auto& d) { return d.count(); });
    }
    
    template <typename Container>
    auto mean_duration(const Container& durations)
    {
        return mean_duration(std::cbegin(durations), std::cend(durations));
    }
    
    template <typename ForwardIt>
    auto stdev_duration(ForwardIt first, ForwardIt last)
    {
        return Maths::stdev(first, last, [] (const auto& d) { return d.count(); });
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
    
    std::string ttc_pad(const std::string& ttc)
    {
        return std::string(16 - ttc.size(), ' ');
    }
    
    template <typename L>
    void print_done(L& log, const TimeInterval& duration)
    {
        const auto time_taken = to_string(duration);
        stream(log) << std::string(15, ' ')
                    << "-"
                    << completed_pad("100%")
                    << "100%"
                    << time_taken_pad(time_taken)
                    << time_taken
                    << ttc_pad("-")
                    << "-";
    }
    
    ProgressMeter::~ProgressMeter()
    {
        if (!done_) {
            const auto now = std::chrono::system_clock::now();
            print_done(log_, TimeInterval {start_, now});
        }
    }
    
    void ProgressMeter::log_completed(const GenomicRegion& completed_region)
    {
        const auto new_bp_processed = right_overhang_size(completed_region, completed_region_);
        
        completed_region_ = encompassing_region(completed_region_, completed_region);
        
        const auto new_percent_done = percent_completed(new_bp_processed, num_bp_to_search_);
        
        percent_unitl_log_ -= new_percent_done;
        
        num_bp_completed_ = region_size(completed_region_);
        
        const auto now = std::chrono::system_clock::now();
        
        if (percent_unitl_log_ <= 0) {
            const auto percent_done = percent_completed(num_bp_completed_, num_bp_to_search_);
            
            const TimeInterval duration {start_, now};
            const auto time_taken = to_string(duration);
            
            if (percent_done >= 100) {
                if (!done_) {
                    print_done(log_, duration);
                    done_ = true;
                }
                return;
            }
            
            const auto percent_since_last_log = percent_done - percent_at_last_log_;
            
            const auto num_blocks_completed = static_cast<std::size_t>(std::floor(percent_since_last_log / percent_block_size_));
            
            const auto duration_since_last_log = std::chrono::duration_cast<DurationUnits>(now - last_log_);
            
            const DurationUnits duration_per_block {static_cast<std::size_t>(duration_since_last_log.count() / num_blocks_completed)};
            
            std::fill_n(std::back_inserter(block_compute_times_),
                        static_cast<std::size_t>(num_blocks_completed),
                        duration_per_block);
            
            const auto num_remaining_blocks = static_cast<std::size_t>((100.0 - percent_done) / percent_block_size_);
            
            remove_outliers(block_compute_times_);
            
            const auto ttc = estimate_ttc(now, block_compute_times_, num_remaining_blocks);
            auto ttc_str   = to_string(ttc);
            
            assert(!ttc_str.empty());
            
            if (ttc_str.front() == '0') {
                ttc_str = "-";
            }
            
            const auto percent_completed = percent_completed_str(num_bp_completed_, num_bp_to_search_);
            
            stream(log_) << curr_position_pad(completed_region)
                         << completed_region.get_contig_name() << ':' << completed_region.get_begin()
                         << completed_pad(percent_completed)
                         << percent_completed
                         << time_taken_pad(time_taken)
                         << time_taken
                         << ttc_pad(ttc_str)
                         << ttc_str;
            
            last_log_            = now;
            percent_unitl_log_   = percent_block_size_;
            percent_at_last_log_ = percent_done;
        }
    }
} // namespace Octopus
