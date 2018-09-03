// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef fixed_mixture_genotype_likelihood_model_hpp
#define fixed_mixture_genotype_likelihood_model_hpp

#include <vector>

#include "core/types/haplotype.hpp"
#include "core/types/genotype.hpp"
#include "core/types/cancer_genotype.hpp"
#include "core/models/haplotype_likelihood_cache.hpp"

namespace octopus { namespace model {

class FixedMixtureGenotypeLikelihoodModel
{
public:
    using LogProbability = double;
    using MixtureVector = std::vector<double>;
    
    FixedMixtureGenotypeLikelihoodModel() = delete;
    
    FixedMixtureGenotypeLikelihoodModel(const HaplotypeLikelihoodCache& likelihoods);
    FixedMixtureGenotypeLikelihoodModel(const HaplotypeLikelihoodCache& likelihoods,
                                        const std::vector<Haplotype>& haplotypes);
    FixedMixtureGenotypeLikelihoodModel(const HaplotypeLikelihoodCache& likelihoods, MixtureVector mixtures);
    FixedMixtureGenotypeLikelihoodModel(const HaplotypeLikelihoodCache& likelihoods, MixtureVector mixtures,
                                        const std::vector<Haplotype>& haplotypes);
    
    FixedMixtureGenotypeLikelihoodModel(const FixedMixtureGenotypeLikelihoodModel&)            = default;
    FixedMixtureGenotypeLikelihoodModel& operator=(const FixedMixtureGenotypeLikelihoodModel&) = default;
    FixedMixtureGenotypeLikelihoodModel(FixedMixtureGenotypeLikelihoodModel&&)                 = default;
    FixedMixtureGenotypeLikelihoodModel& operator=(FixedMixtureGenotypeLikelihoodModel&&)      = default;
    
    ~FixedMixtureGenotypeLikelihoodModel() = default;
    
    const HaplotypeLikelihoodCache& cache() const noexcept;
    const MixtureVector& mixtures() const noexcept;
    
    void prime(const std::vector<Haplotype>& haplotypes);
    void unprime() noexcept;
    bool is_primed() const noexcept;
    
    void set_mixtures(MixtureVector mixtures);
    
    LogProbability evaluate(const Genotype<Haplotype>& genotype) const;
    LogProbability evaluate(const GenotypeIndex& genotype) const;
    LogProbability evaluate(const CancerGenotype<Haplotype>& genotype) const;
    LogProbability evaluate(const CancerGenotypeIndex& genotype) const;

private:
    const HaplotypeLikelihoodCache& likelihoods_;
    MixtureVector mixtures_, log_mixtures_;
    std::vector<HaplotypeLikelihoodCache::LikelihoodVectorRef> indexed_likelihoods_;
    mutable std::vector<HaplotypeLikelihoodCache::LikelihoodVectorRef> likelihood_refs_;
    mutable std::vector<LogProbability> buffer_;
};

template <typename Container1, typename Container2>
Container2&
evaluate(const Container1& genotypes, const FixedMixtureGenotypeLikelihoodModel& model, Container2& result,
         const bool add = false)
{
    if (add) {
        assert(result.size() == genotypes.size());
        std::transform(std::cbegin(genotypes), std::cend(genotypes), std::cbegin(result), std::begin(result),
                       [&] (const auto& genotype, auto curr) { return curr + model.evaluate(genotype); });
    } else {
        result.resize(genotypes.size());
        std::transform(std::cbegin(genotypes), std::cend(genotypes), std::begin(result),
                       [&] (const auto& genotype) { return model.evaluate(genotype); });
    }
    return result;
}

template <typename Container>
auto evaluate(const Container& genotypes, const FixedMixtureGenotypeLikelihoodModel& model)
{
    std::vector<FixedMixtureGenotypeLikelihoodModel::LogProbability> result {};
    evaluate(genotypes, model, result);
    return result;
}

} // namespace model
} // namespace octopus

#endif
