// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef tumour_model_hpp
#define tumour_model_hpp

#include <vector>
#include <unordered_map>
#include <utility>

#include "cancer_genotype_prior_model.hpp"
#include "config/common.hpp"
#include "core/types/haplotype.hpp"
#include "core/models/haplotype_likelihood_cache.hpp"
#include "core/types/cancer_genotype.hpp"

namespace octopus { namespace model {

class TumourModel
{
public:
    struct AlgorithmParameters
    {
        unsigned max_iterations = 1000;
        double epsilon          = 0.05;
    };
    
    struct Priors
    {
        using GenotypeMixturesDirichletAlphas   = std::vector<double>;
        using GenotypeMixturesDirichletAlphaMap = std::unordered_map<SampleName, GenotypeMixturesDirichletAlphas>;
        
        const CancerGenotypePriorModel& genotype_prior_model;
        GenotypeMixturesDirichletAlphaMap alphas;
    };
    
    struct Latents
    {
        using GenotypeMixturesDirichletAlphas   = std::vector<double>;
        using GenotypeMixturesDirichletAlphaMap = std::unordered_map<SampleName, GenotypeMixturesDirichletAlphas>;
        using ProbabilityVector                 = std::vector<double>;
        ProbabilityVector genotype_probabilities;
        GenotypeMixturesDirichletAlphaMap alphas;
    };
    
    struct InferredLatents
    {
        Latents posteriors;
        Latents::ProbabilityVector genotype_log_priors;
        double approx_log_evidence;
    };
    
    TumourModel() = delete;
    
    TumourModel(std::vector<SampleName> samples, Priors priors);
    TumourModel(std::vector<SampleName> samples, Priors priors, AlgorithmParameters parameters);
    
    ~TumourModel() = default;
    
    TumourModel(const TumourModel&)            = default;
    TumourModel& operator=(const TumourModel&) = default;
    TumourModel(TumourModel&&)                 = default;
    TumourModel& operator=(TumourModel&&)      = default;
    
    const Priors& priors() const noexcept;
    
    InferredLatents evaluate(const std::vector<CancerGenotype<Haplotype>>& genotypes,
                             const HaplotypeLikelihoodCache& haplotype_likelihoods) const;
    
    InferredLatents evaluate(const std::vector<CancerGenotype<Haplotype>>& genotypes,
                             const std::vector<CancerGenotypeIndex>& genotype_indices,
                             const HaplotypeLikelihoodCache& haplotype_likelihoods) const;
    
private:
    std::vector<SampleName> samples_;
    Priors priors_;
    AlgorithmParameters parameters_;
};

} // namespace model
} // namespace octopus

#endif
