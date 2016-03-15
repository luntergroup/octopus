//
//  progress_meter.cpp
//  Octopus
//
//  Created by Daniel Cooke on 10/03/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#include "progress_meter.hpp"

#include "mappable_algorithms.hpp"
#include "mappable_map.hpp"
#include "timing.hpp"

namespace Octopus
{
    ProgressMeter::ProgressMeter(const SearchRegions& input_regions)
    :
    num_bp_completed_ {0},
    start_ {std::chrono::system_clock::now()},
    last_log_ {std::chrono::system_clock::now()}
    {
        log_ << "|\t Position \t|\t Complete \t|\t Remaining \t|\t Time taken \t|";
    }
    
    ProgressMeter::ProgressMeter(const GenomicRegion& input_region)
    :
    num_bp_completed_ {0},
    num_bp_until_log_ {0},
    start_ {std::chrono::system_clock::now()},
    last_log_ {std::chrono::system_clock::now()}
    {
        regions_[input_region.get_contig_name()].emplace(input_region);
        
        completed_region_ = head_region(input_region);
        
        num_bp_to_search_ = sum_region_sizes(regions_);
        
        log_ << "|\t Position \t|\t Complete \t|\t Remaining \t|\t Time taken \t|";
    }
    
    double percent_completed(const std::size_t num_bp_completed,
                             const std::size_t num_bp_to_search)
    {
        return 100 * static_cast<double>(num_bp_completed) / num_bp_to_search;
    }
    
    double percent_remaining(const std::size_t num_bp_completed,
                             const std::size_t num_bp_to_search)
    {
        return 100 * (num_bp_to_search - static_cast<double>(num_bp_completed)) / num_bp_to_search;
    }
    
    void ProgressMeter::log_completed(const GenomicRegion& completed_region)
    {
        std::size_t new_bp_processed {0};
        
        if (is_empty_region(completed_region_)) {
            completed_region_ = completed_region;
        } else {
            new_bp_processed = right_overhang_size(completed_region, completed_region_);
            completed_region_ = encompassing_region(completed_region_, completed_region);
        }
        
        num_bp_completed_ = region_size(completed_region_);
        
        last_log_ = std::chrono::system_clock::now();
        
        if (num_bp_until_log_ <= 0) {
            stream(log_) << '\t' << completed_region.get_contig_name() << ':' << completed_region.get_begin()
                                 << "\t\t  "
                                 << percent_completed(num_bp_completed_, num_bp_to_search_) << "%"
                                 << "\t\t  "
                                 << percent_remaining(num_bp_completed_, num_bp_to_search_) << "%"
                                 << "\t\t  "
                                 << TimeInterval {start_, last_log_};
            num_bp_until_log_ = 100000;
        }
        
        num_bp_until_log_ -= new_bp_processed;
    }
} // namespace Octopus
