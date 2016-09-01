// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef read_size_estimator_hpp
#define read_size_estimator_hpp

#include <cstddef>
#include <vector>

#include <boost/optional.hpp>

#include "config/common.hpp"
#include "io/read/read_manager.hpp"

namespace octopus {

boost::optional<std::size_t> estimate_mean_read_size(const std::vector<SampleName>& samples,
                                                     const InputRegionMap& input_regions,
                                                     ReadManager& read_manager,
                                                     unsigned max_sample_size = 1000);

std::size_t default_read_size_estimate() noexcept;

} // namespace octopus

#endif
