//
//  progress_meter.hpp
//  Octopus
//
//  Created by Daniel Cooke on 10/03/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#ifndef progress_meter_hpp
#define progress_meter_hpp

#include <string>
#include <cstddef>
#include <chrono>
#include <deque>

#include "common.hpp"
#include "contig_region.hpp"
#include "genomic_region.hpp"
#include "mappable_flat_set.hpp"
#include "logging.hpp"

namespace Octopus
{
    class ProgressMeter
    {
    public:
        ProgressMeter() = delete;
        
        ProgressMeter(InputRegionMap regions);
        ProgressMeter(GenomicRegion region);
        
        ~ProgressMeter();
        
        ProgressMeter(const ProgressMeter&)            = delete;
        ProgressMeter& operator=(const ProgressMeter&) = delete;
        ProgressMeter(ProgressMeter&&);
        ProgressMeter& operator=(ProgressMeter&&)      = default;
        
        void log_completed(const GenomicRegion& region);
        
    private:
        using RegionSizeType = ContigRegion::SizeType;
        
        using ContigRegionMap = MappableSetMap<ContigNameType, ContigRegion>;
        
        using DurationUnits = std::chrono::milliseconds;
        
        InputRegionMap regions_;
        
        ContigRegionMap completed_regions_;
        
        RegionSizeType num_bp_to_search_, num_bp_completed_;
        
        double percent_block_size_ = 1.0;
        double percent_unitl_log_;
        double percent_at_last_log_;
        
        std::chrono::time_point<std::chrono::system_clock> start_, last_log_;
        
        bool done_;
        
        std::size_t position_tab_length_;
        
        mutable std::deque<DurationUnits> block_compute_times_;
        
        mutable std::mutex mutex_;
        
        Logging::InfoLogger log_;
        
        RegionSizeType merge(const GenomicRegion& region);
        
        void write_header();
        void output_log(const GenomicRegion& region);
        
        std::string curr_position_pad(const GenomicRegion& completed_region) const;
        std::string completed_pad(const std::string& percent_completed) const;
        std::string time_taken_pad(const std::string& time_taken) const;
        std::string ttc_pad(const std::string& ttc) const;
    };
} // namespace Octopus

#endif /* progress_meter_hpp */
