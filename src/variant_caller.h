//
//  variant_caller.h
//  Octopus
//
//  Created by Daniel Cooke on 29/04/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__variant_caller__
#define __Octopus__variant_caller__

#include <vector>
#include <unordered_map>
#include <map>
#include <functional>

#include "common.h"
#include "allele.h"
#include "variant.h"
#include "haplotype_phaser.h"

namespace Octopus
{

namespace VariantCaller
{
    template <typename SampleIdType, typename RealType>
    using AllelePosteriors = std::unordered_map<SampleIdType, std::map<Allele, RealType>>;
    
    AllelePosteriors<Octopus::SampleIdType, Octopus::ProbabilityType>
    get_allele_posteriors(const std::vector<Octopus::SampleIdType>& the_samples,
                          const HaplotypePhaser::PhasedRegions& the_phased_regions,
                          const std::vector<Allele>& the_candidates);
    
    AllelePosteriors<Octopus::SampleIdType, Octopus::ProbabilityType>
    get_allele_posteriors(const std::vector<Octopus::SampleIdType>& the_samples,
                          const HaplotypePhaser::PhasedRegions& the_phased_regions,
                          const std::vector<Variant>& the_candidates);
    
    template <typename SampleIdType, typename RealType>
    using VariantCalls = std::unordered_map<SampleIdType, std::map<Variant, RealType>>;
    
    VariantCalls<Octopus::SampleIdType, Octopus::ProbabilityType>
    call_variants(const std::vector<Octopus::SampleIdType>& the_samples,
                  const HaplotypePhaser::PhasedRegions& the_phased_regions,
                  const std::vector<Variant>& the_candidates);
    
} // end namespace VariantCaller

} // end namespace Octopus

#endif /* defined(__Octopus__variant_caller__) */
