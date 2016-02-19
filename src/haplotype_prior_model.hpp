//
//  haplotype_prior_model.hpp
//  Octopus
//
//  Created by Daniel Cooke on 26/08/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_haplotype_prior_model_hpp
#define Octopus_haplotype_prior_model_hpp

#include <vector>
#include <unordered_map>
#include <functional>
#include <iterator>
#include <cstddef>
#include <cassert>

#include "haplotype.hpp"

namespace Octopus
{
class HaplotypePriorModel
{
public:
    using HaplotypePriorMap = std::unordered_map<std::reference_wrapper<const Haplotype>, double>;
    
    virtual ~HaplotypePriorModel() = default;
    
    double evaluate(const Haplotype& to, const Haplotype& from) const;
    
    HaplotypePriorMap evaluate(std::vector<Haplotype>::const_iterator first,
                               std::vector<Haplotype>::const_iterator last,
                               std::vector<Haplotype>::const_iterator reference);
    
private:
    // ln p(to | from)
    virtual double do_evaluate(const Haplotype& to, const Haplotype& from) const = 0;
    
    virtual HaplotypePriorMap do_evaluate(std::vector<Haplotype>::const_iterator first,
                                          std::vector<Haplotype>::const_iterator last,
                                          std::vector<Haplotype>::const_iterator reference) const = 0;
};

template <typename Container>
HaplotypePriorModel::HaplotypePriorMap
evaluate(const Container& haplotypes, typename Container::const_iterator reference_itr,
         HaplotypePriorModel* prior_model)
{
    assert(reference_itr != std::cend(haplotypes));
    return prior_model->evaluate(std::cbegin(haplotypes), std::cend(haplotypes), reference_itr);
}

void remove_lowest_prior_duplicates(std::vector<Haplotype>& haplotypes,
                                    HaplotypePriorModel::HaplotypePriorMap& haplotype_priors);

namespace debug
{
    void print_haplotype_priors(const HaplotypePriorModel::HaplotypePriorMap& haplotype_priors,
                                std::size_t n = 5);
} // namespace debug

} // namespace Octopus

#endif
