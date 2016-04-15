//
//  cancer_genotype_model.hpp
//  Octopus
//
//  Created by Daniel Cooke on 26/08/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__cancer_genotype_model__
#define __Octopus__cancer_genotype_model__

#include <vector>
#include <unordered_map>
#include <utility>

#include "common.hpp"
#include "haplotype.hpp"
#include "somatic_model.hpp"
#include "haplotype_likelihood_cache.hpp"
#include "cancer_genotype.hpp"

namespace Octopus
{
    namespace GenotypeModel
    {
    class Cancer
    {
    public:
        struct Priors
        {
            using GenotypeMixturesDirichletAlphas   = std::vector<double>;
            using GenotypeMixturesDirichletAlphaMap = std::unordered_map<SampleIdType, GenotypeMixturesDirichletAlphas>;
            
            Priors() = delete;
            template <typename C, typename D> Priors(C&&, D&&);
            ~Priors() = default;
            
            SomaticModel genotype_prior_model;
            GenotypeMixturesDirichletAlphaMap alphas;
        };
        
        struct AlgorithmParameters
        {
            unsigned max_parameter_seeds = 3;
            unsigned max_iterations      = 100;
            double epsilon               = 0.001;
        };
        
        struct Latents
        {
            using GenotypeMixturesDirichletAlphas   = std::vector<double>;
            using GenotypeMixturesDirichletAlphaMap = std::unordered_map<SampleIdType, GenotypeMixturesDirichletAlphas>;
            
            using GenotypeProbabilityMap = std::unordered_map<CancerGenotype<Haplotype>, double>;
            
            Latents() = default;
            template <typename G, typename M> Latents(G&& genotype_probabilities, M&& alphas);
            ~Latents() = default;
            
            GenotypeProbabilityMap genotype_probabilities;
            GenotypeMixturesDirichletAlphaMap alphas;
        };
        
        struct InferredLatents
        {
            InferredLatents(Latents&& posteriors, double approx_log_evidence);
            Latents posteriors;
            double approx_log_evidence;
        };
        
        Cancer() = delete;
        
        explicit Cancer(std::vector<SampleIdType> samples, unsigned ploidy, Priors priors);
        explicit Cancer(std::vector<SampleIdType> samples, unsigned ploidy, Priors priors,
                        AlgorithmParameters parameters);
        
        ~Cancer() = default;
        
        Cancer(const Cancer&)            = default;
        Cancer& operator=(const Cancer&) = default;
        Cancer(Cancer&&)                 = default;
        Cancer& operator=(Cancer&&)      = default;
        
        InferredLatents infer_latents(std::vector<CancerGenotype<Haplotype>> genotypes,
                                      const HaplotypeLikelihoodCache& haplotype_likelihoods) const;
        
    private:
        std::vector<SampleIdType> samples_;
        
        unsigned ploidy_;
        
        Priors priors_;
        
        AlgorithmParameters parameters_;
    };
    
    template <typename C, typename D>
    Cancer::Priors::Priors(C&& genotype_prior_model, D&& alphas)
    :
    genotype_prior_model {std::forward<C>(genotype_prior_model)},
    alphas {std::forward<D>(alphas)}
    {}
    
    template <typename G, typename M>
    Cancer::Latents::Latents(G&& genotype_probabilities, M&& alphas)
    :
    genotype_probabilities {std::forward<G>(genotype_probabilities)},
    alphas {std::forward<M>(alphas)}
    {}
    
    } // namespace GenotypeModel
} // namespace Octopus

#endif /* defined(__Octopus__cancer_genotype_model__) */
