//
//  basic_haplotype_prior_model.cpp
//  Octopus
//
//  Created by Daniel Cooke on 17/02/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#include "basic_haplotype_prior_model.hpp"

#include <algorithm>
#include <iterator>
#include <cassert>

#include "haplotype.hpp"
#include "sequence_utils.hpp"
#include "variant.hpp"
#include "maths.hpp"

#include <iostream> // DEBUG

namespace Octopus
{
BasicHaplotypePriorModel::BasicHaplotypePriorModel(const ReferenceGenome& reference)
:
reference_ {reference}
{}

BasicHaplotypePriorModel::BasicHaplotypePriorModel(const ReferenceGenome& reference,
                                                   double transition_rate, double transversion_rate)
:
reference_ {reference},
transition_rate_ {transition_rate},
transversion_rate_ {transversion_rate}
{}

namespace
{
    std::vector<TandemRepeat> extract_exact_tandem_repeats(const Haplotype& haplotype)
    {
        return find_exact_tandem_repeats(haplotype.get_sequence(), haplotype.get_region(), 1);
    }
} // namespace

double BasicHaplotypePriorModel::do_evaluate(const Haplotype& to, const Haplotype& from) const
{
    const auto mutations = to.difference(from);
    
    double result {1.0};
    
    auto has_indels = std::any_of(std::cbegin(mutations), std::cend(mutations),
                                  [] (const auto& variant) { return is_indel(variant); });
    
    const auto repeats = (has_indels) ? extract_exact_tandem_repeats(from) : std::vector<TandemRepeat> {};
    
    for (const auto& variant : mutations) {
        if (is_snp(variant)) {
            result *= (is_transition(variant)) ? transition_rate_ : transversion_rate_;
        } else if (is_indel(variant)) {
            if (is_insertion(variant)) {
                if (repeats.empty()) {
                    result *= transversion_rate_ * alt_sequence_size(variant);
                } else {
                    result *= transition_rate_ * alt_sequence_size(variant);
                }
            } else {
                if (repeats.empty()) {
                    result *= transversion_rate_ * region_size(variant);
                } else {
                    result *= transition_rate_ * region_size(variant);
                }
            }
        } else {
            auto itr1 = std::cbegin(ref_sequence(variant));
            auto itr2 = std::cbegin(alt_sequence(variant));
            
            std::for_each(itr1, std::cend(ref_sequence(variant)),
                          [this, &itr2, &result] (char base) {
                              if (base != *itr2) {
                                  result *= 0.9 * transversion_rate_;
                              }
                              ++itr2;
                          });
        }
    }
    
    return result;
}

namespace
{
    void normalise(HaplotypePriorModel::HaplotypePriorMap& priors)
    {
        const auto norm = Maths::sum_values(priors);
        for (auto& p : priors) p.second /= norm;
    }
} // namespace

HaplotypePriorModel::HaplotypePriorMap
BasicHaplotypePriorModel::do_compute_maximum_entropy_haplotype_set(std::vector<Haplotype>& haplotypes) const
{
    using std::begin; using std::end; using std::next; using std::distance;
    using std::adjacent_find; using std::for_each; using std::nth_element;
    
    assert(!haplotypes.empty());
    assert(std::is_sorted(std::cbegin(haplotypes), std::cend(haplotypes)));
    
    const Haplotype reference {haplotypes.front().get_region(), reference_.get()};
    
    auto first_duplicate = adjacent_find(begin(haplotypes), end(haplotypes));
    
    if (first_duplicate == end(haplotypes)) {
        HaplotypePriorMap result {haplotypes.size()};
        
        for (const auto& haplotype : haplotypes) {
            result.emplace(haplotype, evaluate(haplotype, reference));
        }
        
        return result;
    }
    
    // need to store copies to prevent iterator invalidation
    std::unordered_map<Haplotype, double, std::hash<Haplotype>, HaveSameAlleles> duplicate_priors {};
    
    duplicate_priors.reserve(distance(first_duplicate, std::end(haplotypes)));
    
    auto it = first_duplicate;
    
    do {
        const auto it2 = std::find_if_not(next(it, 2), end(haplotypes),
                                          [first_duplicate] (const Haplotype& haplotype) {
                                              return haplotype == *first_duplicate;
                                          });
        
        for_each(it, it2, [this, &duplicate_priors, &reference] (const Haplotype& duplicate) {
                duplicate_priors.emplace(duplicate, evaluate(duplicate, reference));
        });
        
        nth_element(it, it, it2, [&duplicate_priors] (const Haplotype& lhs, const Haplotype& rhs) {
                return duplicate_priors.at(lhs) > duplicate_priors.at(rhs);
        });
        
        for_each(next(it), it2, [&duplicate_priors] (const Haplotype& duplicate) {
            return duplicate_priors.erase(duplicate);
        });
        
        it = adjacent_find(it2, end(haplotypes));
    } while (it != end(haplotypes));
    
    haplotypes.erase(std::unique(first_duplicate, end(haplotypes)), end(haplotypes));
    
    HaplotypePriorMap result {haplotypes.size()};
    
    for (const auto& haplotype : haplotypes) {
        if (duplicate_priors.count(haplotype) == 1) {
            result.emplace(haplotype, duplicate_priors.at(haplotype));
        } else {
            result.emplace(haplotype, evaluate(haplotype, reference));
        }
    }
    
    normalise(result);
    
    assert(haplotypes.size() == result.size());
    
    return result;
}

} // namespace Octopus