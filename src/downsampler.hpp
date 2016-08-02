//
//  downsampler.hpp
//  Octopus
//
//  Created by Daniel Cooke on 19/10/2015.
//  Copyright © 2015 Oxford University. All rights reserved.
//

#ifndef downsampler_hpp
#define downsampler_hpp

#include <cstddef>

#include "common.hpp"
#include "mappable_flat_multi_set.hpp"
#include "mappable_map.hpp"
#include "aligned_read.hpp"

namespace octopus { namespace readpipe {

class Downsampler
{
public:
    Downsampler() = default;
    
    Downsampler(unsigned max_coverage, unsigned min_coverage);
    
    Downsampler(const Downsampler&)            = default;
    Downsampler& operator=(const Downsampler&) = default;
    Downsampler(Downsampler&&)                 = default;
    Downsampler& operator=(Downsampler&&)      = default;
    
    ~Downsampler() = default;
    
    // returns the number of reads downsampled
    std::size_t downsample(ReadContainer& reads) const;
    std::size_t downsample(ReadMap& reads) const;
    
private:
    unsigned max_coverage_ = 100'000;
    unsigned min_coverage_ = 100'000;
};

} // namespace readpipe
} // namespace octopus

#endif /* downsampler_hpp */
