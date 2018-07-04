// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef subclone_model_hpp
#define subclone_model_hpp

#include <vector>
#include <unordered_map>
#include <utility>

#include "config/common.hpp"
#include "core/types/haplotype.hpp"
#include "core/types/genotype.hpp"
#include "genotype_prior_model.hpp"
#include "core/models/haplotype_likelihood_cache.hpp"

namespace octopus { namespace model {

class SubcloneModel
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
        
        const GenotypePriorModel& genotype_prior_model;
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
        double approx_log_evidence;
    };
    
    SubcloneModel() = delete;
    
    SubcloneModel(std::vector<SampleName> samples, Priors priors);
    SubcloneModel(std::vector<SampleName> samples, Priors priors, AlgorithmParameters parameters);
    
    SubcloneModel(const SubcloneModel&)            = default;
    SubcloneModel& operator=(const SubcloneModel&) = default;
    SubcloneModel(SubcloneModel&&)                 = default;
    SubcloneModel& operator=(SubcloneModel&&)      = default;
    
    ~SubcloneModel() = default;
    
    const Priors& priors() const noexcept;
    
    InferredLatents evaluate(const std::vector<Genotype<Haplotype>>& genotypes,
                             const HaplotypeLikelihoodCache& haplotype_likelihoods) const;
    
    InferredLatents evaluate(const std::vector<Genotype<Haplotype>>& genotypes,
                             const std::vector<GenotypeIndex>& genotype_indices,
                             const HaplotypeLikelihoodCache& haplotype_likelihoods) const;
    
private:
    std::vector<SampleName> samples_;
    Priors priors_;
    AlgorithmParameters parameters_;
};
    
} // namespace model
} // namespace octopus

#endif
