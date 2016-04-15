//
//  variational_bayes_mixture_model.hpp
//  Octopus
//
//  Created by Daniel Cooke on 12/04/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#ifndef variational_bayes_mixture_model_h
#define variational_bayes_mixture_model_h

#include <array>
#include <vector>

#include "common.hpp"
#include "haplotype.hpp"
#include "haplotype_likelihood_cache.hpp"

namespace Octopus
{
    template <std::size_t K>
    class VariationalBayesMixtureModel
    {
    public:
        struct AlgorithmOptions
        {
            explicit AlgorithmOptions(double epsilon, unsigned max_iterations);
            double epsilon;
            unsigned max_iterations;
        };
        
        struct Latents
        {
            using DirichletAlpha    = std::vector<double>;
            using DirichletAlphas   = std::vector<DirichletAlphas>;
            using ProbabilityVector = std::vector<double>;
            
            template <typename C, typename D> Latents(C&&, D&&);
            
            ProbabilityVector genotype_probabilities;
            DirichletAlphas alphas;
        };
        
        struct InferredLatents
        {
            InferredLatents(Latents&& posteriors, double approx_log_evidence);
            Latents posteriors;
            double approx_log_evidence;
        }
        
        VariationalBayesMixtureModel() = default;
        VariationalBayesMixtureModel(AlgorithmOptions options);
        ~VariationalBayesMixtureModel() = default;
        
        VariationalBayesMixtureModel(const VariationalBayesMixtureModel&)            = default;
        VariationalBayesMixtureModel& operator=(const VariationalBayesMixtureModel&) = default;
        VariationalBayesMixtureModel(VariationalBayesMixtureModel&&)                 = default;
        VariationalBayesMixtureModel& operator=(VariationalBayesMixtureModel&&)      = default;
        
        template <typename Container>
        InferredLatents infer_latents(const std::vector<SampleIdType>& samples,
                                      const Container& genotypes,
                                      const Latents& latents,
                                      const HaplotypeLikelihoodCache& haplotype_log_likelihoods);
        
    private:
        template <std::size_t K>
        using CompressedAlpha = std::array<double, K>;
        
        template <std::size_t K>
        using CompressedAlphas = std::vector<CompressedAlpha<K>>;
        
        class ReadLikelihoods
        {
        public:
            using BaseType = HaplotypeLikelihoodCache::Likelihoods;
            
            ReadLikelihoods() = default;
            explicit ReadLikelihoods(const BaseType&);
            ~ReadLikelihoods() = default;
            
            ReadLikelihoods(const ReadLikelihoods&)            = default;
            ReadLikelihoods& operator=(const ReadLikelihoods&) = default;
            ReadLikelihoods(ReadLikelihoods&&)                 = default;
            ReadLikelihoods& operator=(ReadLikelihoods&&)      = default;
            
            void operator=(const BaseType&);
            void operator=(std::reference_wrapper<const BaseType>);
            
            std::size_t size() const noexcept;
            
            BaseType::const_iterator begin() const noexcept;
            BaseType::const_iterator end() const noexcept;
            
            double operator[](const std::size_t n) const noexcept;
            
        private:
            const BaseType* likelihoods;
        };
        
        template <std::size_t K>
        using CompressedGenotype = std::array<ReadLikelihoods, K>;
        template <std::size_t K>
        using CompressedGenotypes = std::vector<CompressedGenotype<K>>;
        template <std::size_t K>
        using CompressedReadLikelihoods = std::vector<CompressedGenotypes<K>>;
        
        template <std::size_t K>
        using Tau = std::array<double, K>;
        
        template <std::size_t K>
        using ResponsabilityVector = std::vector<Tau<K>>;
        
        template <std::size_t K>
        using ResponsabilityVectors = std::vector<ResponsabilityVector<K>>;
        
        template <std::size_t K>
        struct CompressedLatents
        {
            ProbabilityVector genotype_posteriors;
            LogProbabilityVector genotype_log_posteriors;
            CompressedAlphas<K> alphas;
            ResponsabilityVectors<K> responsabilities;
        };
        
        AlgorithmOptions options_;
    };
} // namespace Octopus

#endif /* variational_bayes_mixture_model_h */
