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

#include "allele.h"
#include "variant.h"

template <typename SampleIdType, typename RealType=double>
using AlleleProbabilityMap = std::unordered_map<SampleIdType, std::map<Allele, RealType>>;

template <typename SampleIdType, typename RealType=double>
using VariantCalls = std::unordered_map<SampleIdType, std::map<Variant, RealType>>;

template <typename SampleIdType, typename RealType=double>
VariantCalls<SampleIdType, RealType>
call_variants(const AlleleProbabilityMap<SampleIdType, RealType>& the_allele_probabilities)
{
    VariantCalls<SampleIdType, RealType> result {};
    
    for (const auto& sample_allele_probabilities : the_allele_probabilities) {
        
    }
    
    return result;
}

#endif /* defined(__Octopus__variant_caller__) */
