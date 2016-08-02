//
//  tumour.hpp
//  Octopus
//
//  Created by Daniel Cooke on 26/08/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__tumour__
#define __Octopus__tumour__

#include <vector>
#include <unordered_map>
#include <utility>

#include "common.hpp"
#include "haplotype.hpp"
#include "somatic_mutation_model.hpp"
#include "haplotype_likelihood_cache.hpp"
#include "cancer_genotype.hpp"

namespace octopus { namespace model
{

class Somatic
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
    
    Somatic() = delete;
    
    Somatic(std::vector<SampleName> samples, unsigned ploidy, Priors priors);
    
    Somatic(std::vector<SampleName> samples, unsigned ploidy, Priors priors,
            AlgorithmParameters parameters);
    
    ~Somatic() = default;
    
    Somatic(const Somatic&)            = default;
    Somatic& operator=(const Somatic&) = default;
    Somatic(Somatic&&)                 = default;
    Somatic& operator=(Somatic&&)      = default;
    
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

#endif /* defined(__Octopus__tumour__) */
