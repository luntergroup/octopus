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

#include "haplotype.hpp"
#include "sequence_utils.hpp"
#include "maths.hpp"
#include "variant.hpp"

#include <iostream> // DEBUG

namespace Octopus
{
BasicHaplotypePriorModel::BasicHaplotypePriorModel(double transition_rate, double transversion_rate)
:
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
    template <typename Map>
    void normalise(Map& haplotype_priors)
    {
        const auto norm = Maths::sum_values(haplotype_priors);
        for (auto& p : haplotype_priors) p.second /= norm;
    }
} // namespace

BasicHaplotypePriorModel::HaplotypePriorMap
BasicHaplotypePriorModel::do_evaluate(std::vector<Haplotype>::const_iterator first,
                                      std::vector<Haplotype>::const_iterator last,
                                      std::vector<Haplotype>::const_iterator reference) const
{
    HaplotypePriorMap result(std::distance(first, last));
    
    std::for_each(first, reference, [&]
                  (const Haplotype& haplotype) {
                      result.emplace(haplotype, evaluate(haplotype, *reference));
                  });
    
    result.emplace(*reference, evaluate(*reference, *reference));
    
    std::for_each(std::next(reference), last,
                  [&] (const Haplotype& haplotype) {
                      result.emplace(haplotype, evaluate(haplotype, *reference));
                  });
    
    normalise(result);
    
    return result;
}
} // namespace Octopus