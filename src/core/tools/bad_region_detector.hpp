// Copyright (c) 2015-2021 Daniel Cooke
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
        enum class Severity { high, low };
        GenomicRegion region;
        Severity severity;
    };
    
    using OptionalReadsReport = boost::optional<const ReadPipe::Report&>;
    
    BadRegionDetector() = default;
    
    BadRegionDetector(Parameters params, boost::optional<const ReadSetProfile&> = boost::none);
    
    BadRegionDetector(const BadRegionDetector&)            = default;
    BadRegionDetector& operator=(const BadRegionDetector&) = default;
    BadRegionDetector(BadRegionDetector&&)                 = default;
    BadRegionDetector& operator=(BadRegionDetector&&)      = default;
    
    ~BadRegionDetector() = default;
    
    std::vector<BadRegion>
    detect(const ReadMap& reads,
           OptionalReadsReport reads_report = boost::none) const;
    
    std::vector<BadRegion>
    detect(const MappableFlatSet<Variant>& variants,
           const ReadMap& reads,
           OptionalReadsReport reads_report = boost::none) const;

private:
    struct InputData
    {
        const ReadMap& reads;
        boost::optional<const MappableFlatSet<Variant>&> variants;
        OptionalReadsReport reads_report;
    };
    
    struct RegionState
    {
        struct ReadSummaryStats
        {
            AlignedRead::NucleotideSequence::size_type max_length;
            std::unordered_map<SampleName, unsigned> average_depths;
            AlignedRead::MappingQuality median_mapping_quality;
        };
        struct VariantSummaryStats
        {
            unsigned count;
            double density;
        };
        GenomicRegion region;
        ReadSummaryStats read_stats;
        boost::optional<VariantSummaryStats> variant_stats;
    };
    
    Parameters params_;
    boost::optional<const ReadSetProfile&> reads_profile_;
    
    std::vector<BadRegion>
    detect(const InputData& data,
           OptionalReadsReport reads_report = boost::none) const;
    std::vector<GenomicRegion>
    find_high_depth_regions(const ReadMap& reads, OptionalReadsReport reads_reports) const;
    std::vector<GenomicRegion>
    get_candidate_bad_regions(const InputData& data) const;
    std::vector<GenomicRegion>
    get_candidate_dense_regions(const MappableFlatSet<Variant>& variants,
                                const ReadMap& reads,
                                OptionalReadsReport reads_report) const;
    double get_max_expected_log_allele_count_per_base() const noexcept;
    void
    fill(RegionState::ReadSummaryStats& stats,
         const ReadMap& reads,
         const GenomicRegion& region,
         OptionalReadsReport reads_report) const;
    void
    fill(RegionState::VariantSummaryStats& stats,
         const MappableFlatSet<Variant>& variants,
         const GenomicRegion& region) const;
    RegionState
    compute_state(const GenomicRegion& region,
                  const InputData& data) const;
    std::vector<RegionState>
    compute_states(const std::vector<GenomicRegion>& region,
                   const InputData& data) const;
    double calculate_probability_good(const RegionState& state, OptionalReadsReport reads_report) const;
    bool is_bad(const RegionState& state, OptionalReadsReport reads_report) const;
};

} // namespace coretools

using coretools::BadRegionDetector;

} // namespace octopus

#endif
