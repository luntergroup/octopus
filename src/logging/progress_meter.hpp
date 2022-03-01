// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef progress_meter_hpp
#define progress_meter_hpp

#include <string>
#include <cstddef>
#include <chrono>
#include <deque>
#include <mutex>

#include "config/common.hpp"
#include "basics/contig_region.hpp"
#include "basics/genomic_region.hpp"
#include "containers/mappable_flat_set.hpp"
#include "logging.hpp"

namespace octopus {

class ProgressMeter
{
public:
    ProgressMeter() = delete;
    
    ProgressMeter(InputRegionMap regions);
    ProgressMeter(GenomicRegion region);
    
    ProgressMeter(const ProgressMeter&)            = delete;
    ProgressMeter& operator=(const ProgressMeter&) = delete;
    ProgressMeter(ProgressMeter&&);
    ProgressMeter& operator=(ProgressMeter&&);
    
    ~ProgressMeter();
    
    void set_max_tick_size(double percent);
    void set_buffer_size(std::size_t n);
    
    void start();
    void resume();
    void pause();
    void stop();
    void reset();
    
    void log_completed(const GenomicRegion& region);
    void log_completed(const GenomicRegion::ContigName& contig);
    
private:
    using RegionSizeType = ContigRegion::Position;
    using ContigRegionMap = MappableSetMap<ContigName, ContigRegion>;
    using DurationUnits = std::chrono::milliseconds;
    
    InputRegionMap target_regions_;
    ContigRegionMap completed_regions_;
    RegionSizeType num_bp_to_search_, num_bp_completed_;
    double min_tick_size_ = 0.1, max_tick_size_ = 1.0, curr_tick_size_ = 1.0;
    double percent_until_tick_;
    double percent_at_last_tick_;
    std::chrono::time_point<std::chrono::system_clock> start_, last_tick_;
    std::deque<std::chrono::duration<float>> tick_durations_;
    bool done_;
    std::size_t position_tab_length_;
    mutable std::deque<DurationUnits> block_compute_times_;
    mutable std::mutex mutex_;
    logging::InfoLogger log_;
    std::size_t buffer_size_;
    std::vector<GenomicRegion> buffer_;
    
    void buffer(const GenomicRegion& region);
    void spill_buffer();
    void do_log_completed(const GenomicRegion& region);
    RegionSizeType merge(const GenomicRegion& region);
    
    void write_header();
    void output_log(const GenomicRegion& region);
    
    std::string position_pad(const GenomicRegion& completed_region) const;
    std::string completed_pad(const std::string& percent_completed, std::size_t position_tick_size) const;
    std::string time_taken_pad(const std::string& time_taken) const;
    std::string ttc_pad(const std::string& ttc) const;
    void update_tick_size();
};

} // namespace octopus

#endif
