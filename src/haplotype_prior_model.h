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
#include <algorithm> // std::transform, std::for_each
#include <numeric>   // std::accumulate

#include "variant.h"
#include "haplotype.h"

template <typename T>
double get_haplotype_prior_probability(const Haplotype& the_haplotype,
                                       const std::vector<Variant>& possible_variants)
{
    T result {1};
    
    for (const auto& variant : possible_variants) {
        if (the_haplotype.contains(variant)) result *= variant.get_segregation_probability();
    }
    
    return result;
}

template <typename T>
std::vector<T> get_haplotype_prior_pseudo_counts(const std::vector<Haplotype>& the_haplotypes,
                                                 const std::vector<Variant>& possible_variants,
                                                 T beta, T c)
{
    std::vector<T> result(the_haplotypes.size());
    
    std::transform(the_haplotypes.cbegin(), the_haplotypes.cend(), result.begin(),
                   [&possible_variants] (const auto& haplotype) {
                       return get_haplotype_prior_probability(haplotype, possible_variants);
                   });
    
    auto total = std::accumulate(result.cbegin(), result.cend(), 0);
    
    std::for_each(result.begin(), result.end(), [total, beta, c] (auto p) { beta * p /= total + c; });
    
    return result;
}

#endif /* defined(__Octopus__haplotype_prior_model__) */
