// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef input_reads_profiler_hpp
#define input_reads_profiler_hpp

#include <cstddef>
#include <vector>

#include <boost/optional.hpp>

#include "config/common.hpp"
#include "io/read/read_manager.hpp"
#include "readpipe/read_pipe.hpp"

namespace octopus {

struct ReadSetProfileConfig
{
    unsigned max_samples_per_sample = 10;
    unsigned max_sample_size = 1000;
};

struct ReadSetProfile
{
    std::size_t mean_read_bytes, read_bytes_stdev;
    std::size_t mean_depth, depth_stdev;
    std::vector<std::size_t> sample_mean_depth;
    std::vector<std::size_t> sample_depth_stdev;
};

boost::optional<ReadSetProfile> profile_reads(const std::vector<SampleName>& samples,
                                              const InputRegionMap& input_regions,
                                              const ReadManager& source,
                                              ReadSetProfileConfig config = ReadSetProfileConfig {});

boost::optional<std::size_t> estimate_mean_read_size(const std::vector<SampleName>& samples,
                                                     const InputRegionMap& input_regions,
                                                     ReadManager& read_manager,
                                                     unsigned max_sample_size = 1000);

std::size_t default_read_size_estimate() noexcept;

} // namespace octopus

#endif
