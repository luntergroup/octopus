//
//  dirichlet_model.cpp
//  Octopus
//
//  Created by Daniel Cooke on 17/02/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#include "dirichlet_model.hpp"

#include "map_utils.hpp"

namespace Octopus
{
namespace GenotypeModel
{
    HaplotypeFrequencyMap init_haplotype_frequencies(const std::vector<Haplotype>& haplotypes)
    {
        HaplotypeFrequencyMap result {haplotypes.size()};
        
        const double uniform {1.0 / haplotypes.size()};
        
        for (const auto& haplotype : haplotypes) {
            result.emplace(haplotype, uniform);
        }
        
        return result;
    }
    
    HaplotypeFrequencyMap init_haplotype_frequencies(const HaplotypePriorCountMap& haplotype_counts)
    {
        HaplotypeFrequencyMap result {haplotype_counts.size()};
        
        auto n = Maths::sum_values(haplotype_counts);
        
        for (const auto& haplotype_count : haplotype_counts) {
            result.emplace(haplotype_count.first, haplotype_count.second / n);
        }
        
        return result;
    }
    
    HaplotypePriorCountMap
    compute_haplotype_prior_counts(const HaplotypeFrequencyMap& haplotype_priors)
    {
        static constexpr double   PRECISION      {40.0};
        static constexpr unsigned MAX_ITERATIONS {100};
        
        HaplotypePriorCountMap result {};
        
        if (haplotype_priors.empty()) return result;
        
        const auto alphas = Maths::dirichlet_mle(extract_values(haplotype_priors),
                                                 PRECISION, MAX_ITERATIONS);
        
        result.reserve(haplotype_priors.size());
        
        std::transform(std::cbegin(haplotype_priors), std::cend(haplotype_priors),
                       std::cbegin(alphas), std::inserter(result, begin(result)),
                       [] (const auto& p, const double alpha) {
                           return std::make_pair(std::ref(p.first), alpha);
                       });
        
        for (auto& h : result) h.second = 1; // DEBUG - uniform priors
        
        return result;
    }
} // namespace GenotypeModel
} // namespace Octopus