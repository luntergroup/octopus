// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef tumour_model_hpp
#define tumour_model_hpp

#include <vector>
#include <unordered_map>
#include <utility>

#include "config/common.hpp"
#include "core/types/haplotype.hpp"
#include "core/models/haplotype_likelihood_cache.hpp"
#include "core/types/cancer_genotype.hpp"
#include "core/models/mutation/somatic_mutation_model.hpp"

namespace octopus { namespace model
{

class TumourModel
{
public:
    struct AlgorithmParameters
    {
        unsigned max_parameter_seeds = 3;
        unsigned max_iterations      = 100;
        double epsilon               = 0.001;
    };
    
    struct Priors
    {
        using GenotypeMixturesDirichletAlphas   = std::vector<double>;
        using GenotypeMixturesDirichletAlphaMap = std::unordered_map<SampleName, GenotypeMixturesDirichletAlphas>;
        
        SomaticMutationModel genotype_prior_model;
        GenotypeMixturesDirichletAlphaMap alphas;
    };
    
    struct Latents
    {
        using GenotypeMixturesDirichletAlphas   = std::vector<double>;
        using GenotypeMixturesDirichletAlphaMap = std::unordered_map<SampleName, GenotypeMixturesDirichletAlphas>;
        
        using GenotypeProbabilityMap = std::unordered_map<CancerGenotype<Haplotype>, double>;
        
        GenotypeProbabilityMap genotype_probabilities;
        GenotypeMixturesDirichletAlphaMap alphas;
    };
    
    struct InferredLatents
    {
        Latents posteriors;
        double approx_log_evidence;
    };
    
    TumourModel() = delete;
    
    TumourModel(std::vector<SampleName> samples, unsigned ploidy, Priors priors);
    
    TumourModel(std::vector<SampleName> samples, unsigned ploidy, Priors priors,
                AlgorithmParameters parameters);
    
    ~TumourModel() = default;
    
    TumourModel(const TumourModel&)            = default;
    TumourModel& operator=(const TumourModel&) = default;
    TumourModel(TumourModel&&)                 = default;
    TumourModel& operator=(TumourModel&&)      = default;
    
    InferredLatents infer_latents(std::vector<CancerGenotype<Haplotype>> genotypes,
                                  const HaplotypeLikelihoodCache& haplotype_likelihoods) const;
    
private:
    std::vector<SampleName> samples_;
    unsigned ploidy_;
    Priors priors_;
    AlgorithmParameters parameters_;
};

} // namespace model
} // namespace octopus

#endif
