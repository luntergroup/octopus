//
//  cnv.hpp
//  Octopus
//
//  Created by Daniel Cooke on 12/04/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#ifndef cnv_hpp
#define cnv_hpp

#include <vector>
#include <unordered_map>
#include <utility>

#include "common.hpp"
#include "haplotype.hpp"
#include "genotype.hpp"
#include "coalescent_model.hpp"
#include "haplotype_likelihood_cache.hpp"

namespace octopus { namespace model {

class CNV
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
    
    CNV() = delete;
    
    CNV(std::vector<SampleName> samples, unsigned ploidy, Priors priors);
    
    CNV(std::vector<SampleName> samples, unsigned ploidy, Priors priors,
        AlgorithmParameters parameters);
    
    CNV(const CNV&)            = default;
    CNV& operator=(const CNV&) = default;
    CNV(CNV&&)                 = default;
    CNV& operator=(CNV&&)      = default;
    
    ~CNV() = default;
    
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

#endif /* cnv_hpp */
