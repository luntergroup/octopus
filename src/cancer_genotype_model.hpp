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
#include "coalescent_model.hpp"
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
            
            template <typename C, typename D> Priors(C&&, D&&);
            
            CoalescentModel genotype_prior_model;
            GenotypeMixturesDirichletAlphaMap alphas;
        };
        
        struct AlgorithmParameters
        {
            unsigned max_parameter_seeds = 3;
            unsigned max_iterations = 200;
            double epsilon = 0.01;
        };
        
        struct InferredLatents
        {
            using GenotypeMixturesDirichletAlphas   = std::vector<double>;
            using GenotypeMixturesDirichletAlphaMap = std::unordered_map<SampleIdType, GenotypeMixturesDirichletAlphas>;
            using GenotypePosteriorMap = std::unordered_map<CancerGenotype<Haplotype>, double>;
            
            InferredLatents() = default;
            template <typename G, typename M> InferredLatents(G&& genotype_posteriors, M&& alphas);
            
            GenotypePosteriorMap genotype_posteriors;
            GenotypeMixturesDirichletAlphaMap alphas;
        };
        
        explicit Cancer(std::vector<SampleIdType> samples, const SampleIdType& normal_sample,
                        unsigned ploidy, Priors priors);
        explicit Cancer(std::vector<SampleIdType> samples, const SampleIdType& normal_sample,
                        unsigned ploidy, Priors priors, AlgorithmParameters parameters);
        
        InferredLatents infer_latents(std::vector<CancerGenotype<Haplotype>> genotypes,
                                      const HaplotypeLikelihoodCache& haplotype_likelihoods) const;
        
        InferredLatents infer_latents(const std::vector<Haplotype>& haplotypes,
                                      const HaplotypeLikelihoodCache& haplotype_likelihoods) const;
        
    private:
        std::vector<SampleIdType> samples_;
        SampleIdType normal_sample_;
        
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
    Cancer::InferredLatents::InferredLatents(G&& genotype_posteriors, M&& alphas)
    :
    genotype_posteriors {std::forward<G>(genotype_posteriors)},
    alphas {std::forward<M>(alphas)}
    {}
    
    } // namespace GenotypeModel
} // namespace Octopus

#endif /* defined(__Octopus__cancer_genotype_model__) */
