// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef variable_mixture_genotype_likelihood_model_hpp
#define variable_mixture_genotype_likelihood_model_hpp

#include <vector>

#include "core/types/haplotype.hpp"
#include "core/types/genotype.hpp"
#include "core/types/cancer_genotype.hpp"
#include "core/models/haplotype_likelihood_array.hpp"

namespace octopus { namespace model {

class VariableMixtureGenotypeLikelihoodModel
{
public:
    using LogProbability = double;
    using MixtureVector = std::vector<double>;
    
    VariableMixtureGenotypeLikelihoodModel() = delete;
    
    VariableMixtureGenotypeLikelihoodModel(const HaplotypeLikelihoodArray& likelihoods);
    VariableMixtureGenotypeLikelihoodModel(const HaplotypeLikelihoodArray& likelihoods,
                                           const std::vector<Haplotype>& haplotypes);
    VariableMixtureGenotypeLikelihoodModel(const HaplotypeLikelihoodArray& likelihoods, MixtureVector mixtures);
    VariableMixtureGenotypeLikelihoodModel(const HaplotypeLikelihoodArray& likelihoods, MixtureVector mixtures,
                                           const std::vector<Haplotype>& haplotypes);
    
    VariableMixtureGenotypeLikelihoodModel(const VariableMixtureGenotypeLikelihoodModel&)            = default;
    VariableMixtureGenotypeLikelihoodModel& operator=(const VariableMixtureGenotypeLikelihoodModel&) = default;
    VariableMixtureGenotypeLikelihoodModel(VariableMixtureGenotypeLikelihoodModel&&)                 = default;
    VariableMixtureGenotypeLikelihoodModel& operator=(VariableMixtureGenotypeLikelihoodModel&&)      = default;
    
    ~VariableMixtureGenotypeLikelihoodModel() = default;
    
    const HaplotypeLikelihoodArray& cache() const noexcept;
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
    const HaplotypeLikelihoodArray& likelihoods_;
    MixtureVector mixtures_, log_mixtures_;
    std::vector<HaplotypeLikelihoodArray::LikelihoodVectorRef> indexed_likelihoods_;
    mutable std::vector<HaplotypeLikelihoodArray::LikelihoodVectorRef> likelihood_refs_;
    mutable std::vector<HaplotypeLikelihoodArray::LogProbability> buffer_;
};

template <typename Container1, typename Container2>
Container2&
evaluate(const Container1& genotypes, const VariableMixtureGenotypeLikelihoodModel& model, Container2& result,
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
auto evaluate(const Container& genotypes, const VariableMixtureGenotypeLikelihoodModel& model)
{
    std::vector<VariableMixtureGenotypeLikelihoodModel::LogProbability> result {};
    evaluate(genotypes, model, result);
    return result;
}

} // namespace model
} // namespace octopus

#endif
