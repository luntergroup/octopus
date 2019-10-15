// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef input_reads_profiler_hpp
#define input_reads_profiler_hpp

#include <cstddef>
#include <vector>
#include <iosfwd>

#include <boost/optional.hpp>

#include "config/common.hpp"
#include "basics/aligned_read.hpp"
#include "io/read/read_manager.hpp"
#include "readpipe/read_pipe.hpp"

namespace octopus {

struct ReadSetProfileConfig
{
    unsigned max_draws_per_sample = 500;
    unsigned max_reads_per_draw = 10'000;
    boost::optional<AlignedRead::NucleotideSequence::size_type> fragment_size = boost::none;
};

struct ReadSetProfile
{
    template <typename T>
    struct SummaryStats
    {
        T max, min, mean, median, stdev;
    };
    
    using DepthSummaryStats = SummaryStats<std::size_t>;
    struct DepthStats
    {
        using DiscreteDistribution = std::vector<double>;
        DiscreteDistribution distribution;
        DepthSummaryStats all, positive;
    };
    struct GenomeContigDepthStatsPair
    {
        using ContigDepthStatsMap = std::unordered_map<GenomicRegion::ContigName, DepthStats>;
        ContigDepthStatsMap contig;
        DepthStats genome;
    };
    struct SampleCombinedDepthStatsPair
    {
        std::unordered_map<SampleName, GenomeContigDepthStatsPair> sample;
        GenomeContigDepthStatsPair combined;
    };
    using MappingQualityStats = SummaryStats<AlignedRead::MappingQuality>;
    using ReadLengthStats = SummaryStats<AlignedRead::NucleotideSequence::size_type>;
    using ReadMemoryStats = SummaryStats<MemoryFootprint>;
    
    SampleCombinedDepthStatsPair depth_stats;
    ReadMemoryStats memory_stats;
    boost::optional<ReadMemoryStats> fragmented_memory_stats;
    ReadLengthStats length_stats;
    MappingQualityStats mapping_quality_stats;
};

boost::optional<ReadSetProfile>
profile_reads(const std::vector<SampleName>& samples,
              const InputRegionMap& input_regions,
              const ReadManager& source,
              ReadSetProfileConfig config = ReadSetProfileConfig {});

std::ostream& operator<<(std::ostream& os, const ReadSetProfile& profile);

} // namespace octopus

#endif
