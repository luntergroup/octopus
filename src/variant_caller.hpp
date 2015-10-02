//
//  variant_caller.hpp
//  Octopus
//
//  Created by Daniel Cooke on 26/08/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_variant_caller_hpp
#define Octopus_variant_caller_hpp

#include <vector>
#include <iterator>
#include <algorithm> // std::max
#include <limits>    // std::numeric_limits
#include <cmath>     // std::log10

#include "common.hpp"
#include "reference_genome.hpp"
#include "read_manager.hpp"
#include "read_filter.hpp"
#include "read_transform.hpp"
#include "candidate_variant_generator.hpp"

class GenomicRegion;
class Variant;
class VcfRecord;

namespace Octopus
{
    class VariantCaller
    {
    public:
        using ReadFilter = ReadFilter<Octopus::ReadContainer::const_iterator>;
        
        VariantCaller() = delete;
        VariantCaller(ReferenceGenome& reference, ReadManager& read_manager,
                      ReadFilter read_filter, ReadTransform read_transform,
                      CandidateVariantGenerator& candidate_generator);
        virtual ~VariantCaller() = default;
        
        VariantCaller(const VariantCaller&)            = delete;
        VariantCaller& operator=(const VariantCaller&) = delete;
        VariantCaller(VariantCaller&&)                 = delete;
        VariantCaller& operator=(VariantCaller&&)      = delete;
        
        std::vector<VcfRecord> call_variants(const GenomicRegion& region);
        
    protected:
        using ReadMap = Octopus::ReadMap;
        
        ReferenceGenome& reference_;
        ReadManager& read_manager_;
        ReadFilter read_filter_;
        ReadTransform read_transform_;
        CandidateVariantGenerator& candidate_generator_;
        
    private:
        bool done_calling(const GenomicRegion& region) const noexcept;
        
        virtual GenomicRegion get_init_region(const GenomicRegion& region) = 0;
        virtual GenomicRegion get_next_region(const GenomicRegion& current_region) = 0;
        virtual std::vector<VcfRecord> call_variants(const GenomicRegion& region,
                                                     const std::vector<Variant>& candidates,
                                                     const ReadMap& reads) = 0;
    };
    
    inline unsigned to_phred_quality(double p)
    {
        return static_cast<unsigned>(-10.0 * std::log10(std::max(1.0 - p, std::numeric_limits<double>::epsilon())));
    }
    
} // namespace Octopus

#endif
