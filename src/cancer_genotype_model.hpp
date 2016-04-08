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
#include <array>
#include <utility>

#include "common.hpp"
#include "haplotype.hpp"
#include "haplotype_likelihood_cache.hpp"
#include "cancer_genotype.hpp"
#include "probability_matrix.hpp"

namespace Octopus
{
    namespace GenotypeModel
    {
    class Cancer
    {
    public:
        static constexpr unsigned K {3};
        
        using GenotypeMixturesDirichletAlphas   = std::array<double, K>;
        using GenotypeMixturesDirichletAlphaMap = std::unordered_map<SampleIdType, GenotypeMixturesDirichletAlphas>;
        using GenotypePosteriorMap              = ProbabilityMatrix<Genotype<Haplotype>>;
        
        struct Priors
        {
            GenotypeMixturesDirichletAlphaMap alphas;
        };
        
        struct AlgorithmParameters
        {
            unsigned max_parameter_seeds = 3;
            unsigned max_iterations = 100;
            double epsilon = 0.001;
        };
        
        struct Latents
        {
            Latents() = default;
            template <typename G, typename M>
            Latents(G&& genotype_posteriors, M&& alphas);
            
            GenotypePosteriorMap genotype_posteriors;
            GenotypeMixturesDirichletAlphaMap alphas;
        };
        
        explicit Cancer(std::vector<SampleIdType> samples, const SampleIdType& normal_sample,
                        Priors priors);
        explicit Cancer(std::vector<SampleIdType> samples, const SampleIdType& normal_sample,
                        Priors priors, AlgorithmParameters parameters);
        
        Latents infer_latents(const std::vector<Haplotype>& haplotypes,
                              const HaplotypeLikelihoodCache& haplotype_likelihoods) const;
        
    private:
        std::vector<SampleIdType> samples_;
        SampleIdType normal_sample_;
        
        Priors priors_;
        
        AlgorithmParameters parameters_;
    };
    
    template <typename G, typename M>
    Cancer::Latents::Latents(G&& genotype_posteriors, M&& alphas)
    :
    genotype_posteriors {std::forward<G>(genotype_posteriors)},
    alphas {std::forward<M>(alphas)}
    {}
    
    } // namespace GenotypeModel
} // namespace Octopus

#endif /* defined(__Octopus__cancer_genotype_model__) */
