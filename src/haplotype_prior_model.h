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
#include <algorithm> // std::count_if
#include <unordered_map>
#include <boost/math/distributions/poisson.hpp>

#include "allele.h"
#include "variant.h"
#include "haplotype.h"
#include "maths.h"

#include <iostream> // TEST

namespace Octopus
{

template <typename RealType, typename ForwardIterator>
RealType prior_probability(const Haplotype& the_haplotype, ForwardIterator first_possible_variant,
                           ForwardIterator last_possible_variant)
{
    auto num_non_reference_alleles = std::count_if(first_possible_variant, last_possible_variant,
                                                   [&the_haplotype] (const Variant& variant) {
                                                       return the_haplotype.contains(variant.get_alternative_allele());
                                                   });
    
    boost::math::poisson_distribution<RealType> poisson {0.003 * size(the_haplotype)};
    
    //return std::pow(0.0003333, num_non_reference_alleles);
    
    return boost::math::pdf(poisson, num_non_reference_alleles) / 3;
}

template <typename RealType, typename Haplotypes, typename ForwardIterator>
std::unordered_map<Haplotype, RealType>
get_haplotype_prior_probabilities(const Haplotypes& the_haplotypes, ForwardIterator first_possible_variant,
                                  ForwardIterator last_possible_variant)
{
    std::unordered_map<Haplotype, RealType> result {};
    result.reserve(the_haplotypes.size());
    
    for (const auto& haplotype : the_haplotypes) {
        result.emplace(haplotype, prior_probability<RealType>(haplotype, first_possible_variant,
                                                              last_possible_variant));
    }
    
    return result;
}

} // end namespace Octopus

#endif /* defined(__Octopus__haplotype_prior_model__) */
