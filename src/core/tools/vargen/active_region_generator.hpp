// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef active_region_generator_hpp
#define active_region_generator_hpp

#include <string>

#include "config/common.hpp"
#include "basics/genomic_region.hpp"
#include "basics/aligned_read.hpp"
#include "io/reference/reference_genome.hpp"

namespace octopus { namespace coretools {

class ActiveRegionGenerator
{
public:
    struct Options {};
    
    ActiveRegionGenerator() = delete;
    
    ActiveRegionGenerator(const ReferenceGenome& reference);
    ActiveRegionGenerator(const ReferenceGenome& reference, Options options);
    
    ActiveRegionGenerator(const ActiveRegionGenerator&)            = default;
    ActiveRegionGenerator& operator=(const ActiveRegionGenerator&) = default;
    ActiveRegionGenerator(ActiveRegionGenerator&&)                 = default;
    ActiveRegionGenerator& operator=(ActiveRegionGenerator&&)      = default;
    
    ~ActiveRegionGenerator() = default;
    
    void add_read(const SampleName& sample, const AlignedRead& read);
    template <typename ForwardIterator>
    void add_reads(const SampleName& sample, ForwardIterator first_read, ForwardIterator last_read);
    
    std::vector<GenomicRegion> generate(const GenomicRegion& region, const std::string& generator) const;
    
    void clear() noexcept;
    
private:

};

template <typename ForwardIterator>
void ActiveRegionGenerator::add_reads(const SampleName& sample, ForwardIterator first_read, ForwardIterator last_read)
{

}

} // namespace coretools
} // namespace octopus

#endif
