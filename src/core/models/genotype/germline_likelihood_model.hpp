// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef germline_likelihood_model_hpp
#define germline_likelihood_model_hpp

#include <vector>

#include "core/types/haplotype.hpp"
#include "core/types/genotype.hpp"
#include "core/models/haplotype_likelihood_array.hpp"

namespace octopus { namespace model {

class GermlineLikelihoodModel
{
public:
    using LogProbability = double;
    
    GermlineLikelihoodModel() = delete;
    
    GermlineLikelihoodModel(const HaplotypeLikelihoodArray& likelihoods);
    GermlineLikelihoodModel(const HaplotypeLikelihoodArray& likelihoods, const std::vector<Haplotype>& haplotypes);
    
    GermlineLikelihoodModel(const GermlineLikelihoodModel&)            = default;
    GermlineLikelihoodModel& operator=(const GermlineLikelihoodModel&) = default;
    GermlineLikelihoodModel(GermlineLikelihoodModel&&)                 = default;
    GermlineLikelihoodModel& operator=(GermlineLikelihoodModel&&)      = default;
    
    ~GermlineLikelihoodModel() = default;
    
    const HaplotypeLikelihoodArray& cache() const noexcept;
    
    void prime(const std::vector<Haplotype>& haplotypes);
    void unprime() noexcept;
    bool is_primed() const noexcept;
    
    LogProbability evaluate(const Genotype<Haplotype>& genotype) const;
    LogProbability evaluate(const GenotypeIndex& genotype) const;
    
private:
    const HaplotypeLikelihoodArray& likelihoods_;
    std::vector<HaplotypeLikelihoodArray::LikelihoodVectorRef> indexed_likelihoods_;
    mutable std::vector<HaplotypeLikelihoodArray::LikelihoodVectorRef> likelihood_refs_;
    mutable std::vector<LogProbability> buffer_;
    
    // These are just for optimisation
    LogProbability evaluate_haploid(const Genotype<Haplotype>& genotype) const;
    LogProbability evaluate_diploid(const Genotype<Haplotype>& genotype) const;
    LogProbability evaluate_triploid(const Genotype<Haplotype>& genotype) const;
    LogProbability evaluate_tetraploid(const Genotype<Haplotype>& genotype) const;
    LogProbability evaluate_polyploid(const Genotype<Haplotype>& genotype) const;
};

template <typename Container1, typename Container2>
Container2&
evaluate(const Container1& genotypes, const GermlineLikelihoodModel& model, Container2& result)
{
    result.resize(genotypes.size());
    std::transform(std::cbegin(genotypes), std::cend(genotypes), std::begin(result),
                   [&] (const auto& genotype) { return model.evaluate(genotype); });
    return result;
}

template <typename Container>
auto evaluate(const Container& genotypes, const GermlineLikelihoodModel& model)
{
    std::vector<GermlineLikelihoodModel::LogProbability> result {};
    evaluate(genotypes, model, result);
    return result;
}

} // namespace model
} // namespace octopus

#endif
