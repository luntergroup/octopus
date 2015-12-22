//
//  downsampler.hpp
//  Octopus
//
//  Created by Daniel Cooke on 19/10/2015.
//  Copyright © 2015 Oxford University. All rights reserved.
//

#ifndef downsampler_hpp
#define downsampler_hpp

#include "common.hpp"
#include "mappable_set.hpp"
#include "mappable_map.hpp"
#include "aligned_read.hpp"

namespace Octopus
{

MappableSet<AlignedRead>
downsample(const MappableSet<AlignedRead>& reads, unsigned max_coverage, unsigned min_coverage);

template <typename T>
MappableMap<T, AlignedRead>
downsample(const MappableMap<T, AlignedRead>& reads, const unsigned max_coverage, const unsigned min_coverage)
{
    MappableMap<T, AlignedRead> result {};
    result.reserve(reads.size());
    
    for (const auto& sample_reads : reads) {
        result.emplace(sample_reads.first, downsample(sample_reads.second, max_coverage, min_coverage));
    }
    
    return result;
}

class Downsampler
{
public:
    Downsampler() = default;
    Downsampler(unsigned max_coverage, unsigned min_coverage);
    ~Downsampler() = default;
    
    Downsampler(const Downsampler&)            = default;
    Downsampler& operator=(const Downsampler&) = default;
    Downsampler(Downsampler&&)                 = default;
    Downsampler& operator=(Downsampler&&)      = default;
    
    template <typename R>
    R operator()(R&& reads)
    {
        return downsample(std::forward<R>(reads), max_coverage_, min_coverage_);
    }
    
private:
    unsigned max_coverage_ = 100'000;
    unsigned min_coverage_ = 100'000;
};

} // namespace Octopus

#endif /* downsampler_hpp */
