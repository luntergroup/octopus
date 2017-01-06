// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef variational_bayes_mixture_model_hpp
#define variational_bayes_mixture_model_hpp

#include <array>
#include <vector>
#include <functional>

#include "config/common.hpp"
#include "core/types/haplotype.hpp"
#include "core/models/haplotype_likelihood_cache.hpp"

namespace octopus { namespace model {

template <std::size_t K>
class VariationalBayesMixtureModel
{
public:
    using DirichletAlpha  = std::vector<double>;
    using DirichletAlphas = std::vector<DirichletAlphas>;
    
    template <typename Model>
    struct Priors
    {
        Model genotype_model;
        DirichletAlphas alphas;
    };
    
    struct Posteriors
    {
        using ProbabilityVector = std::vector<double>;
        ProbabilityVector genotype_probabilities;
        DirichletAlphas alphas;
        double log_evidence;
    };
    
    struct Options
    {
        double epsilon;
        unsigned max_iterations;
    };
    
    using Seed = std::vector<double>;
    using SeedGenerator = std::function<Seed()>;
    using SeedGeneratorSet = std::vector<SeedGenerator>;
    
    VariationalBayesMixtureModel() = default;
    
    VariationalBayesMixtureModel(SeedGeneratorSet seed_generators, Options options);
    
    VariationalBayesMixtureModel(const VariationalBayesMixtureModel&)            = default;
    VariationalBayesMixtureModel& operator=(const VariationalBayesMixtureModel&) = default;
    VariationalBayesMixtureModel(VariationalBayesMixtureModel&&)                 = default;
    VariationalBayesMixtureModel& operator=(VariationalBayesMixtureModel&&)      = default;
    
    ~VariationalBayesMixtureModel() = default;
    
    template <typename G>
    Posteriors infer_latents(const std::vector<SampleName>& samples,
                             const std::vector<G>& genotypes,
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
        
        ReadLikelihoods(const ReadLikelihoods&)            = default;
        ReadLikelihoods& operator=(const ReadLikelihoods&) = default;
        ReadLikelihoods(ReadLikelihoods&&)                 = default;
        ReadLikelihoods& operator=(ReadLikelihoods&&)      = default;
    
        ~ReadLikelihoods() = default;
        
        void operator=(const BaseType&);
        void operator=(std::reference_wrapper<const BaseType>);
        
        std::size_t size() const noexcept;
        
        BaseType::const_iterator begin() const noexcept;
        BaseType::const_iterator end() const noexcept;
        
        double operator[](const std::size_t n) const noexcept;
    
    private:
        const BaseType *likelihoods;
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
    
    SeedGeneratorSet seed_generators_;
    Options options_;
};
    
} // namespace model
} // namespace octopus

#endif
