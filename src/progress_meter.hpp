//
//  progress_meter.hpp
//  Octopus
//
//  Created by Daniel Cooke on 10/03/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#ifndef progress_meter_hpp
#define progress_meter_hpp

#include <cstddef>
#include <chrono>
#include <deque>

#include "common.hpp"
#include "genomic_region.hpp"
#include "logging.hpp"

namespace Octopus
{
    class ProgressMeter
    {
    public:
        ProgressMeter() = delete;
        
        ProgressMeter(const SearchRegions& input_regions);
        ProgressMeter(const GenomicRegion& input_region);
        
        ~ProgressMeter();
        
        ProgressMeter(const ProgressMeter&)            = delete;
        ProgressMeter& operator=(const ProgressMeter&) = delete;
        ProgressMeter(ProgressMeter&&)                 = delete;
        ProgressMeter& operator=(ProgressMeter&&)      = delete;
        
        void log_completed(const GenomicRegion& completed_region);
        
    private:
        SearchRegions regions_;
        
        std::size_t num_bp_to_search_, num_bp_completed_;
        
        GenomicRegion completed_region_;
        
        double percent_block_size_ = 1.0;
        double percent_unitl_log_;
        double percent_at_last_log_;
        
        std::chrono::time_point<std::chrono::system_clock> start_, last_log_;
        
        bool done_;
        
        Logging::InfoLogger log_;
        
        using DurationUnits = std::chrono::milliseconds;
        
        std::deque<DurationUnits> block_compute_times_;
        
        mutable std::mutex mutex_;
    };
} // namespace Octopus

#endif /* progress_meter_hpp */
