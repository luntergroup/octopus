//
//  haplotype_prior_model.cpp
//  Octopus
//
//  Created by Daniel Cooke on 01/10/2015.
//  Copyright © 2015 Oxford University. All rights reserved.
//

#include "haplotype_prior_model.hpp"

#include <algorithm>
#include <iostream>

namespace Octopus
{
double HaplotypePriorModel::evaluate(const Haplotype& to, const Haplotype& from) const
{
    return this->do_evaluate(to, from);
}

HaplotypePriorModel::HaplotypePriorMap
HaplotypePriorModel::evaluate(std::vector<Haplotype>::const_iterator first,
                              std::vector<Haplotype>::const_iterator last,
                              std::vector<Haplotype>::const_iterator reference)
{
    return this->do_evaluate(first, last, reference);
}

// non-member methods

void remove_lowest_prior_duplicates(std::vector<Haplotype>& haplotypes,
                                    HaplotypePriorModel::HaplotypePriorMap& haplotype_priors)
{
    assert(std::is_sorted(std::cbegin(haplotypes), std::cend(haplotypes)));
    
    auto first_duplicate = std::begin(haplotypes);
    auto last_duplicate  = first_duplicate;
    
    while (true) {
        first_duplicate = std::adjacent_find(first_duplicate, std::end(haplotypes));
        
        if (first_duplicate == std::end(haplotypes)) break;
        
        last_duplicate = std::find_if_not(std::next(first_duplicate, 2), std::end(haplotypes),
                                          [=] (const Haplotype& haplotype) {
                                              return haplotype == *first_duplicate;
                                          });
        
        std::nth_element(first_duplicate, first_duplicate, last_duplicate,
                         [&] (const Haplotype& lhs, const Haplotype& rhs) {
                             return haplotype_priors.at(lhs) > haplotype_priors.at(rhs);
                         });
        
        std::cout << "removing " << std::distance(std::next(first_duplicate), last_duplicate) << " duplicates" << std::endl;
        
        std::for_each(std::next(first_duplicate), last_duplicate,
                      [&] (const Haplotype& haplotype) {
                          haplotype_priors.erase(haplotype);
                      });
        
        first_duplicate = last_duplicate;
    }
    
    haplotype_priors.rehash(haplotype_priors.size());
    
    haplotypes.erase(std::unique(std::begin(haplotypes), std::end(haplotypes)), std::end(haplotypes));
    
    assert(haplotypes.size() == haplotype_priors.size());
}

namespace debug
{
    void print_haplotype_priors(const HaplotypePriorModel::HaplotypePriorMap& haplotype_priors,
                                const std::size_t n)
    {
        auto m = std::min(haplotype_priors.size(), n);
        
        std::cout << "printing top " << m << " haplotype priors" << std::endl;
        
        using HaplotypeReference = std::reference_wrapper<const Haplotype>;
        
        std::vector<std::pair<HaplotypeReference, double>> v {};
        v.reserve(haplotype_priors.size());
        
        std::copy(std::cbegin(haplotype_priors), std::cend(haplotype_priors), std::back_inserter(v));
        
        std::sort(std::begin(v), std::end(v),
                  [] (const auto& lhs, const auto& rhs) {
                      return lhs.second > rhs.second;
                  });
        
        for (unsigned i {0}; i < m; ++i) {
            print_variant_alleles(v[i].first);
            std::cout << " " << std::setprecision(10) << v[i].second << std::endl;
        }
    }
} // namespace debug
} // namespace Octopus
