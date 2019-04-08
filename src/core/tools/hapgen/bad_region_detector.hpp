// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef bad_region_detector_hpp
#define bad_region_detector_hpp

#include <vector>

#include <boost/optional.hpp>

#include "config/common.hpp"
#include "basics/genomic_region.hpp"
#include "core/types/variant.hpp"
#include "readpipe/read_pipe.hpp"
#include "utils/input_reads_profiler.hpp"

namespace octopus { namespace coretools {

class BadRegionDetector
{
public:
    struct Parameters
    {
        enum class Tolerance { low, normal, high };
        double heterozygosity, heterozygosity_stdev;
        Tolerance tolerance = Tolerance::normal;
    };
    
    struct BadRegion
    {
        enum class RecommendedAction { skip, restrict_lagging };
        GenomicRegion region;
        RecommendedAction action;
    };
    
    BadRegionDetector() = default;
    
    BadRegionDetector(Parameters params, boost::optional<ReadSetProfile> = boost::none);
    
    BadRegionDetector(const BadRegionDetector&)            = default;
    BadRegionDetector& operator=(const BadRegionDetector&) = default;
    BadRegionDetector(BadRegionDetector&&)                 = default;
    BadRegionDetector& operator=(BadRegionDetector&&)      = default;
    
    ~BadRegionDetector() = default;
    
    std::vector<BadRegion>
    detect(const MappableFlatSet<Variant>& candidate_variants, const ReadMap& reads,
           boost::optional<const ReadPipe::Report&> reads_report = boost::none) const;

private:
    struct RegionState
    {
        GenomicRegion region;
        unsigned variant_count;
        double variant_density;
        AlignedRead::NucleotideSequence::size_type max_read_length;
        std::unordered_map<SampleName, unsigned> sample_mean_read_depths;
        AlignedRead::MappingQuality median_mapping_quality;
    };
    
    Parameters params_;
    boost::optional<ReadSetProfile> reads_profile_;
    
    std::vector<GenomicRegion>
    get_candidate_bad_regions(const MappableFlatSet<Variant>& candidate_variants, const ReadMap& reads,
                              boost::optional<const ReadPipe::Report&> reads_report) const;
    std::vector<GenomicRegion>
    get_candidate_dense_regions(const MappableFlatSet<Variant>& candidate_variants, const ReadMap& reads,
                                boost::optional<const ReadPipe::Report&> reads_report) const;
    double get_max_expected_log_allele_count_per_base() const noexcept;
    RegionState
    compute_state(const GenomicRegion& region, const MappableFlatSet<Variant>& candidate_variants,
                  const ReadMap& reads, boost::optional<const ReadPipe::Report&> reads_report) const;
    std::vector<RegionState>
    compute_states(const std::vector<GenomicRegion>& region, const MappableFlatSet<Variant>& candidate_variants,
                   const ReadMap& reads, boost::optional<const ReadPipe::Report&> reads_report) const;
    double calculate_probability_good(const RegionState& state) const;
    bool is_bad(const RegionState& state) const;
};

} // namespace coretools
} // namespace octopus

#endif
