//
//  haplotype_prior_model.cpp
//  Octopus
//
//  Created by Daniel Cooke on 01/10/2015.
//  Copyright © 2015 Oxford University. All rights reserved.
//

#include <haplotype_prior_model.hpp>

#include <algorithm> // std::any_of

#include "haplotype.hpp"
#include "variant_utils.hpp"
#include "sequence_utils.hpp"
#include "maths.hpp"

namespace Octopus
{

HaplotypePriorModel::HaplotypePriorModel(double transition_rate, double transversion_rate)
:
transition_rate_ {transition_rate},
transversion_rate_ {transversion_rate}
{}

std::vector<TandemRepeat> find_exact_tandem_repeats(const Haplotype& haplotype)
{
    return find_exact_tandem_repeats(haplotype.get_sequence(), haplotype.get_region(), 1);
}

double HaplotypePriorModel::evaluate(const Haplotype& to, const Haplotype& from)
{
    auto mutations = to.difference(from);
    
    double result {1.0};
    
    auto has_indels = std::any_of(std::cbegin(mutations), std::cend(mutations),
                                  [] (const auto& variant) { return is_indel(variant); });
    
    auto repeats = (has_indels) ? find_exact_tandem_repeats(from) : std::vector<TandemRepeat> {};
    
    for (const auto& variant : mutations) {
        if (is_snp(variant)) {
            result *= (is_transition(variant)) ? transition_rate_ : transversion_rate_;
        } else if (is_indel(variant)) {
            if (is_insertion(variant)) {
                if (repeats.empty()) {
                    result *= transversion_rate_ * variant.alternative_allele_size();
                } else {
                    result *= transition_rate_ * variant.alternative_allele_size();
                }
            } else {
                if (repeats.empty()) {
                    result *= transversion_rate_ * size(variant);
                } else {
                    result *= transition_rate_ * size(variant);
                }
            }
        } else {
            auto itr1 = std::cbegin(variant.get_reference_allele_sequence());
            auto itr2 = std::cbegin(variant.get_alternative_allele_sequence());
            
            std::for_each(itr1, std::cend(variant.get_reference_allele_sequence()),
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

std::unordered_map<Haplotype, double>
HaplotypePriorModel::evaluate(const std::vector<Haplotype>& haplotypes, const Haplotype& reference)
{
    std::unordered_map<Haplotype, double> result {};
    result.reserve(haplotypes.size());
    
    for (const auto& haplotype : haplotypes) {
        result.emplace(haplotype, evaluate(haplotype, reference));
    }
    
    const auto norm = Maths::sum_values(result);
    
    for (auto& p : result) p.second /= norm;
    
    return result;
}

} // namespace Octopus
