//
//  cnv_model.hpp
//  octopus
//
//  Created by Daniel Cooke on 12/04/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#ifndef cnv_model_hpp
#define cnv_model_hpp

#include <vector>
#include <unordered_map>
#include <utility>

#include <config/common.hpp>
#include <core/types/haplotype.hpp>
#include <core/types/genotype.hpp>
#include <core/models/genotype/coalescent_model.hpp>
#include <core/models/haplotype_likelihood_cache.hpp>

namespace octopus { namespace model {

class CNVModel
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
        
        CoalescentModel genotype_prior_model;
        GenotypeMixturesDirichletAlphaMap alphas;
    };
    
    struct Latents
    {
        using GenotypeMixturesDirichletAlphas   = std::vector<double>;
        using GenotypeMixturesDirichletAlphaMap = std::unordered_map<SampleName, GenotypeMixturesDirichletAlphas>;
        
        using GenotypeProbabilityMap = std::unordered_map<Genotype<Haplotype>, double>;
        
        GenotypeProbabilityMap genotype_probabilities;
        GenotypeMixturesDirichletAlphaMap alphas;
    };
    
    struct InferredLatents
    {
        Latents posteriors;
        double approx_log_evidence;
    };
    
    CNVModel() = delete;
    
    CNVModel(std::vector<SampleName> samples, unsigned ploidy, Priors priors);
    
    CNVModel(std::vector<SampleName> samples, unsigned ploidy, Priors priors,
             AlgorithmParameters parameters);
    
    CNVModel(const CNVModel&)            = default;
    CNVModel& operator=(const CNVModel&) = default;
    CNVModel(CNVModel&&)                 = default;
    CNVModel& operator=(CNVModel&&)      = default;
    
    ~CNVModel() = default;
    
    InferredLatents infer_latents(std::vector<Genotype<Haplotype>> genotypes,
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
