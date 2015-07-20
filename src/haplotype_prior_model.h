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
#include "variant_utils.h"

#include <iostream> // TEST

namespace Octopus
{

namespace detail
{
    // A ts/tv ratio of around 2.1 seems to be a good approximation for human genome.
    // See https://www.biostars.org/p/4751/ for a discussion on this.
    template <typename RealType> constexpr RealType transition_rate {0.000222}; //0.0006999 or 0.000222?
    template <typename RealType> constexpr RealType transversion_rate {0.000111}; //0.0003333 or 0.000111?
}

template <typename RealType, typename ForwardIterator>
RealType prior_probability(const Haplotype& the_haplotype, ForwardIterator first_possible_variant,
                           ForwardIterator last_possible_variant)
{
    auto num_transitions = std::count_if(first_possible_variant, last_possible_variant,
                                         [&the_haplotype] (const Variant& variant) {
                                             return is_transition(variant) && the_haplotype.contains(variant.get_alternative_allele());
                                         });
    
    auto num_transversions = std::count_if(first_possible_variant, last_possible_variant,
                                           [&the_haplotype] (const Variant& variant) {
                                               return is_transversion(variant) && the_haplotype.contains(variant.get_alternative_allele());
                                           });
    
//    auto num_insertions = std::count_if(first_possible_variant, last_possible_variant,
//                                        [&the_haplotype] (const Variant& variant) {
//                                            return is_insertion(variant) && the_haplotype.contains(variant.get_alternative_allele());
//                                        });
//    
//    auto num_deletions = std::count_if(first_possible_variant, last_possible_variant,
//                                       [&the_haplotype] (const Variant& variant) {
//                                           return is_deletion(variant) && the_haplotype.contains(variant.get_alternative_allele());
//                                       });
    
//    auto num_non_reference_alleles = num_transitions + num_transversions + num_insertions + num_deletions;
//    
//    boost::math::poisson_distribution<RealType> transiton_poisson {
//        detail::transition_rate<RealType> * size(the_haplotype)
//    };
//    
//    boost::math::poisson_distribution<RealType> transversion_poisson {
//        detail::transversion_rate<RealType> * size(the_haplotype)
//    };
    
    constexpr RealType r {detail::transition_rate<RealType> + detail::transversion_rate<RealType>};
    
    return std::pow(1 - r, size(the_haplotype) - (num_transitions + num_transversions)) *
            std::pow(detail::transition_rate<RealType>, num_transitions) *
            std::pow(detail::transversion_rate<RealType>, num_transversions);
    
    //return boost::math::pdf(transiton_poisson, num_transitions + num_transversions);
}

template <typename RealType, typename Container, typename ForwardIterator>
std::unordered_map<Haplotype, RealType>
get_haplotype_prior_probabilities(const Container& haplotypes, ForwardIterator first_possible_variant,
                                  ForwardIterator last_possible_variant)
{
    std::unordered_map<Haplotype, RealType> result {};
    result.reserve(haplotypes.size());
    
    for (const auto& haplotype : haplotypes) {
        result.emplace(haplotype, prior_probability<RealType>(haplotype, first_possible_variant,
                                                              last_possible_variant));
    }
    
    auto norm = sum_values(result);
    
    for (auto& p : result) {
        p.second /= norm;
    }
    
    return result;
}

} // end namespace Octopus

#endif /* defined(__Octopus__haplotype_prior_model__) */
