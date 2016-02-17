//
//  dirichlet_model.cpp
//  Octopus
//
//  Created by Daniel Cooke on 17/02/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#include "dirichlet_model.hpp"

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
    
    HaplotypePriorCountMap compute_haplotype_prior_counts(const std::vector<Haplotype>& haplotypes,
                                                          const ReferenceGenome& reference,
                                                          const HaplotypePriorModel& haplotype_prior_model)
    {
        using std::begin; using std::cbegin; using std::cend; using std::transform;
        
        HaplotypePriorCountMap result {};
        
        if (haplotypes.empty()) return result;
        
        const Haplotype reference_haplotype {haplotypes.front().get_region(), reference};
        
        std::vector<double> p(haplotypes.size());
        transform(cbegin(haplotypes), cend(haplotypes), begin(p),
                  [&haplotype_prior_model, &reference_haplotype] (const auto& haplotype) {
                      return haplotype_prior_model.evaluate(haplotype, reference_haplotype);
                  });
        
        const auto norm = std::accumulate(cbegin(p), cend(p), 0.0);
        
        transform(cbegin(p), cend(p), begin(p), [norm] (auto x) { return x / norm; });
        
        constexpr double precision {40.0};
        constexpr unsigned max_iterations {100};
        const auto alphas = Maths::dirichlet_mle(p, precision, max_iterations);
        
        result.reserve(haplotypes.size());
        
        transform(cbegin(haplotypes), cend(haplotypes), cbegin(alphas), std::inserter(result, begin(result)),
                  [] (const auto& haplotype, double a) { return std::make_pair(std::ref(haplotype), a); });
        
        for (auto& h : result) h.second = 1; // DEBUG - uniform priors
        
        return result;
    }
} // namespace GenotypeModel
} // namespace Octopus