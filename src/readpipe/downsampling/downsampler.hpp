// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef downsampler_hpp
#define downsampler_hpp

#include <cstddef>

#include "config/common.hpp"
#include "containers/mappable_flat_multi_set.hpp"
#include "containers/mappable_map.hpp"
#include "basics/aligned_read.hpp"

namespace octopus { namespace readpipe {

/**
 Downsampler removes reads from a collection in order to reduce the coverage at positions
 which exceed a maximum bound.
 */
class Downsampler
{
public:
    Downsampler() = default;
    
    Downsampler(unsigned trigger_coverage, unsigned target_coverage);
    
    Downsampler(const Downsampler&)            = default;
    Downsampler& operator=(const Downsampler&) = default;
    Downsampler(Downsampler&&)                 = default;
    Downsampler& operator=(Downsampler&&)      = default;
    
    ~Downsampler() = default;
    
    // Returns the number of reads removed
    std::size_t downsample(ReadContainer& reads) const;
    
private:
    unsigned trigger_coverage_ = 10'000;
    unsigned target_coverage_  = 10'000;
};

template <typename Map>
auto downsample(Map& reads, const Downsampler& downsampler)
{
    std::size_t result {0};
    
    for (auto& p : reads) {
        result += downsampler.downsample(p.second);
    }
    
    return result;
}

} // namespace readpipe
} // namespace octopus

#endif
