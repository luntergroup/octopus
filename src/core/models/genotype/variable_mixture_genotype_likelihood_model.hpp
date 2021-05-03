// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef variable_mixture_genotype_likelihood_model_hpp
#define variable_mixture_genotype_likelihood_model_hpp

#include <vector>

#include "core/types/haplotype.hpp"
#include "core/types/indexed_haplotype.hpp"
#include "core/types/genotype.hpp"
#include "core/types/cancer_genotype.hpp"
#include "core/types/partitioned_genotype.hpp"
#include "core/models/haplotype_likelihood_array.hpp"

namespace octopus { namespace model {

class VariableMixtureGenotypeLikelihoodModel
{
public:
    using LogProbability = double;
    using MixtureVector = std::vector<double>;
    
    VariableMixtureGenotypeLikelihoodModel() = delete;
    
    VariableMixtureGenotypeLikelihoodModel(const HaplotypeLikelihoodArray& likelihoods);
    VariableMixtureGenotypeLikelihoodModel(const HaplotypeLikelihoodArray& likelihoods, MixtureVector mixtures);
    
    VariableMixtureGenotypeLikelihoodModel(const VariableMixtureGenotypeLikelihoodModel&)            = default;
    VariableMixtureGenotypeLikelihoodModel& operator=(const VariableMixtureGenotypeLikelihoodModel&) = delete;
    VariableMixtureGenotypeLikelihoodModel(VariableMixtureGenotypeLikelihoodModel&&)                 = default;
    VariableMixtureGenotypeLikelihoodModel& operator=(VariableMixtureGenotypeLikelihoodModel&&)      = delete;
    
    ~VariableMixtureGenotypeLikelihoodModel() = default;
    
    const HaplotypeLikelihoodArray& cache() const noexcept;
    const MixtureVector& mixtures() const noexcept;
    
    void set_mixtures(MixtureVector mixtures);
    
    LogProbability evaluate(const Genotype<Haplotype>& genotype) const;
    LogProbability evaluate(const Genotype<IndexedHaplotype<>>& genotype) const;
    LogProbability evaluate(const CancerGenotype<Haplotype>& genotype) const;
    LogProbability evaluate(const CancerGenotype<IndexedHaplotype<>>& genotype) const;
    template <typename IndexType, std::size_t K>
    LogProbability evaluate(const PartitionedGenotype<IndexedHaplotype<IndexType>, K>& genotype) const;

private:
    const HaplotypeLikelihoodArray& likelihoods_;
    MixtureVector mixtures_, log_mixtures_;
    mutable std::vector<HaplotypeLikelihoodArray::LikelihoodVectorRef> likelihood_refs_;
    mutable std::vector<HaplotypeLikelihoodArray::LogProbability> buffer_;
};

template <typename IndexType, std::size_t K>
VariableMixtureGenotypeLikelihoodModel::LogProbability 
VariableMixtureGenotypeLikelihoodModel::evaluate(const PartitionedGenotype<IndexedHaplotype<IndexType>, K>& genotype) const
{
    assert(genotype.ploidy() == mixtures_.size());
    assert(buffer_.size() == mixtures_.size());
    LogProbability result {0};
    const auto num_likelihoods = likelihoods_.num_likelihoods();
    for (std::size_t read_idx {0}; read_idx < num_likelihoods; ++read_idx) {
        const auto get_likelihood = [=] (const auto& haplotype, auto log_mixture) noexcept {
                           return log_mixture + likelihoods_[haplotype][read_idx]; };
        auto buffer_itr = std::begin(buffer_);
        unsigned offset {0};
        for (std::size_t k {0}; k < K; ++k) {
            const auto& partition = genotype.partition(k);
            buffer_itr = std::transform(std::cbegin(partition), std::cend(partition),
                                        std::next(std::cbegin(log_mixtures_), offset), 
                                        buffer_itr,
                                        get_likelihood);
            offset += partition.ploidy();
        }
        result += maths::log_sum_exp(buffer_);
    }
    return result;
}

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
