//
//  haplotype_prior_model.h
//  Octopus
//
//  Created by Daniel Cooke on 06/04/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__haplotype_prior_model__
#define __Octopus__haplotype_prior_model__

#include <vector>

#include "allele.h"
#include "variant.h"
#include "haplotype.h"
#include "maths.h"

template <typename RealType>
RealType get_haplotype_prior_probability(const Haplotype& the_haplotype,
                                         const std::vector<Variant>& possible_variants)
{
    RealType result {1};
    
    for (const auto& variant : possible_variants) {
        if (the_haplotype.contains(variant)) result *= variant.get_segregation_probability();
    }
    
    return result;
}

#endif /* defined(__Octopus__haplotype_prior_model__) */
