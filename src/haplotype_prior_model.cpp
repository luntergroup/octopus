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
