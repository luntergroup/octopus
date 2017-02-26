// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef assembler_active_region_generator_hpp
#define assembler_active_region_generator_hpp

#include <vector>
#include <functional>

#include "utils/coverage_tracker.hpp"
#include "io/reference/reference_genome.hpp"

namespace octopus {

class GenomicRegion;
class AlignedRead;

namespace coretools {

class AssemblerActiveRegionGenerator
{
public:
    enum class TriggerType { snv, indel, structual };
    
    AssemblerActiveRegionGenerator() = delete;
    
    AssemblerActiveRegionGenerator(const ReferenceGenome& reference);
    AssemblerActiveRegionGenerator(const ReferenceGenome& reference, std::vector<TriggerType> trigger_types);
    
    AssemblerActiveRegionGenerator(const AssemblerActiveRegionGenerator&)            = default;
    AssemblerActiveRegionGenerator& operator=(const AssemblerActiveRegionGenerator&) = default;
    AssemblerActiveRegionGenerator(AssemblerActiveRegionGenerator&&)                 = default;
    AssemblerActiveRegionGenerator& operator=(AssemblerActiveRegionGenerator&&)      = default;
    
    ~AssemblerActiveRegionGenerator() = default;
    
    void add(const AlignedRead& read);
    
    std::vector<GenomicRegion> generate(const GenomicRegion& region) const;

private:
    std::reference_wrapper<const ReferenceGenome> reference_;
    bool snvs_interesting_ = false, indels_interesting_ = true, structual_interesting_ = false;
    CoverageTracker coverage_tracker_, interesting_read_coverages_, clipped_coverage_tracker_;
    
    bool is_interesting(const AlignedRead& read) const;
};
    
} // namespace coretools
} // namespace octopus

#endif

