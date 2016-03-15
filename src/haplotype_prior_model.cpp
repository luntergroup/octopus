//
//  haplotype_prior_model.cpp
//  Octopus
//
//  Created by Daniel Cooke on 01/10/2015.
//  Copyright © 2015 Oxford University. All rights reserved.
//

#include "haplotype_prior_model.hpp"

#include <iostream>

#include "maths.hpp"

namespace Octopus
{
double HaplotypePriorModel::evaluate(const Haplotype& to, const Haplotype& from) const
{
    return this->do_evaluate(to, from);
}

HaplotypePriorModel::HaplotypePriorMap
HaplotypePriorModel::compute_maximum_entropy_haplotype_set(std::vector<Haplotype>& haplotypes) const
{
    return this->do_compute_maximum_entropy_haplotype_set(haplotypes);
}

// non-member methods

namespace
{
    template <typename Map>
    void normalise(Map& haplotype_priors)
    {
        const auto norm = Maths::sum_values(haplotype_priors);
        for (auto& p : haplotype_priors) p.second /= norm;
    }
} // namespace

namespace debug
{
    void print_haplotype_priors(const HaplotypePriorModel::HaplotypePriorMap& haplotype_priors,
                                const std::size_t n)
    {
        print_haplotype_priors(std::cout, haplotype_priors, n);
    }
} // namespace debug
} // namespace Octopus
