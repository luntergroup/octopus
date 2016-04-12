//
//  cnv_genotype_model.hpp
//  Octopus
//
//  Created by Daniel Cooke on 12/04/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#ifndef cnv_genotype_model_hpp
#define cnv_genotype_model_hpp

#include <vector>
#include <unordered_map>
#include <utility>

#include "common.hpp"
#include "haplotype.hpp"
#include "genotype.hpp"
#include "coalescent_model.hpp"
#include "haplotype_likelihood_cache.hpp"

namespace Octopus
{
namespace GenotypeModel
{
    class CNV
    {
    public:
        struct Priors
        {
            using GenotypeMixturesDirichletAlphas   = std::vector<double>;
            using GenotypeMixturesDirichletAlphaMap = std::unordered_map<SampleIdType, GenotypeMixturesDirichletAlphas>;
            
            Priors() = delete;
            template <typename C, typename D> Priors(C&&, D&&);
            ~Priors() = default;
            
            CoalescentModel genotype_prior_model;
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
            
            using GenotypeProbabilityMap = std::unordered_map<Genotype<Haplotype>, double>;
            
            Latents() = default;
            template <typename G, typename M> Latents(G&& genotype_posteriors, M&& alphas);
            ~Latents() = default;
            
            GenotypeProbabilityMap genotype_posteriors;
            GenotypeMixturesDirichletAlphaMap alphas;
        };
        
        struct InferredLatents
        {
            InferredLatents(Latents&& posterior_latents, double approx_log_evidence);
            Latents posterior_latents;
            double approx_log_evidence;
        };
        
        CNV() = delete;
        explicit CNV(std::vector<SampleIdType> samples, unsigned ploidy, Priors priors);
        explicit CNV(std::vector<SampleIdType> samples, unsigned ploidy, Priors priors,
                     AlgorithmParameters parameters);
        ~CNV() = default;
        
        CNV(const CNV&)            = default;
        CNV& operator=(const CNV&) = default;
        CNV(CNV&&)                 = default;
        CNV& operator=(CNV&&)      = default;
        
        InferredLatents infer_latents(std::vector<Genotype<Haplotype>> genotypes,
                                      const HaplotypeLikelihoodCache& haplotype_likelihoods) const;
        
    private:
        std::vector<SampleIdType> samples_;
        
        unsigned ploidy_;
        
        Priors priors_;
        
        AlgorithmParameters parameters_;
    };
    
    template <typename C, typename D>
    CNV::Priors::Priors(C&& genotype_prior_model, D&& alphas)
    :
    genotype_prior_model {std::forward<C>(genotype_prior_model)},
    alphas {std::forward<D>(alphas)}
    {}
    
    template <typename G, typename M>
    CNV::Latents::Latents(G&& genotype_posteriors, M&& alphas)
    :
    genotype_posteriors {std::forward<G>(genotype_posteriors)},
    alphas {std::forward<M>(alphas)}
    {}
    
} // namespace GenotypeModel
} // namespace Octopus

#endif /* cnv_genotype_model_hpp */
