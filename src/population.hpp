//
//  population.hpp
//  Octopus
//
//  Created by Daniel Cooke on 26/08/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__population__
#define __Octopus__population__

#include <vector>
#include <unordered_map>
#include <functional>

#include "common.hpp"
#include "haplotype.hpp"
#include "coalescent_model.hpp"
#include "haplotype_likelihood_cache.hpp"
#include "genotype.hpp"
#include "probability_matrix.hpp"
#include "logging.hpp"

namespace octopus { namespace model {

class Population
{
public:
    struct Latents
    {
        using GenotypeProbabilityVector       = std::vector<double>;
        using SampleGenotypeProbabilityVector = std::vector<GenotypeProbabilityVector>;
        
        SampleGenotypeProbabilityVector genotype_probabilities;
    };
    
    struct InferredLatents
    {
        using HaplotypeReference    = std::reference_wrapper<const Haplotype>;
        using HaplotypePosteriorMap = std::unordered_map<HaplotypeReference, double>;
        
        Latents posteriors;
        HaplotypePosteriorMap haplotype_posteriors;
        double log_evidence;
    };
    
    using GenotypeVector = std::vector<Genotype<Haplotype>>;
    
    Population() = delete;
    
    Population(const CoalescentModel& genotype_prior_model);
    
    ~Population() = default;
    
    Population(const Population&)            = delete;
    Population& operator=(const Population&) = delete;
    Population(Population&&)                 = delete;
    Population& operator=(Population&&)      = delete;
    
    InferredLatents infer_latents(const std::vector<SampleName>& samples,
                                  const GenotypeVector& genotypes,
                                  const std::vector<Haplotype>& haplotypes,
                                  const HaplotypeLikelihoodCache& haplotype_likelihoods) const;
    
private:
    std::reference_wrapper<const CoalescentModel> genotype_prior_model_;
};

} // namesapce model
} // namespace octopus

#endif /* defined(__Octopus__population__) */
