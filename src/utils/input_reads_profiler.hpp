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
    unsigned max_samples_per_sample = 10;
    unsigned max_sample_size = 10'000;
    boost::optional<AlignedRead::NucleotideSequence::size_type> fragment_size = boost::none;
};

struct ReadSetProfile
{
    std::size_t mean_read_bytes, read_bytes_stdev;
    std::size_t mean_depth, median_depth, depth_stdev, median_positive_depth, mean_positive_depth;
    std::vector<SampleName> samples;
    std::vector<std::size_t> sample_mean_depth, sample_median_depth, sample_median_positive_depth, sample_mean_positive_depth;
    std::vector<std::size_t> sample_depth_stdev;
    AlignedRead::NucleotideSequence::size_type max_read_length, median_read_length;
    AlignedRead::MappingQuality max_mapping_quality, median_mapping_quality, rmq_mapping_quality;
    boost::optional<std::size_t> fragmented_template_median_bytes;
};

boost::optional<ReadSetProfile>
profile_reads(const std::vector<SampleName>& samples,
              const InputRegionMap& input_regions,
              const ReadManager& source,
              ReadSetProfileConfig config = ReadSetProfileConfig {});

std::ostream& operator<<(std::ostream& os, const ReadSetProfile& profile);

} // namespace octopus

#endif
