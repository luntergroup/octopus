// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef population_model_hpp
#define population_model_hpp

#include <vector>
#include <unordered_map>
#include <functional>

#include "config/common.hpp"
#include "core/types/haplotype.hpp"
#include "core/models/genotype/coalescent_model.hpp"
#include "core/models/haplotype_likelihood_cache.hpp"
#include "core/types/genotype.hpp"
#include "containers/probability_matrix.hpp"
#include "logging/logging.hpp"

namespace octopus { namespace model {

class PopulationModel
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
    
    PopulationModel() = delete;
    
    PopulationModel(const CoalescentModel& genotype_prior_model);
    
    PopulationModel(const PopulationModel&)            = delete;
    PopulationModel& operator=(const PopulationModel&) = delete;
    PopulationModel(PopulationModel&&)                 = delete;
    PopulationModel& operator=(PopulationModel&&)      = delete;
    
    ~PopulationModel() = default;
    
    InferredLatents infer_latents(const std::vector<SampleName>& samples,
                                  const GenotypeVector& genotypes,
                                  const std::vector<Haplotype>& haplotypes,
                                  const HaplotypeLikelihoodCache& haplotype_likelihoods) const;
    
private:
    std::reference_wrapper<const CoalescentModel> genotype_prior_model_;
};

} // namesapce model
} // namespace octopus

#endif
