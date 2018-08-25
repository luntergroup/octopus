// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef germline_likelihood_model_hpp
#define germline_likelihood_model_hpp

#include <vector>

#include "core/types/haplotype.hpp"
#include "core/types/genotype.hpp"
#include "core/models/haplotype_likelihood_cache.hpp"

namespace octopus { namespace model {

class GermlineLikelihoodModel
{
public:
    GermlineLikelihoodModel() = delete;
    
    GermlineLikelihoodModel(const HaplotypeLikelihoodCache& likelihoods);
    GermlineLikelihoodModel(const HaplotypeLikelihoodCache& likelihoods, const std::vector<Haplotype>& haplotypes);
    
    GermlineLikelihoodModel(const GermlineLikelihoodModel&)            = default;
    GermlineLikelihoodModel& operator=(const GermlineLikelihoodModel&) = default;
    GermlineLikelihoodModel(GermlineLikelihoodModel&&)                 = default;
    GermlineLikelihoodModel& operator=(GermlineLikelihoodModel&&)      = default;
    
    ~GermlineLikelihoodModel() = default;
    
    const HaplotypeLikelihoodCache& cache() const noexcept;
    
    void prime(const std::vector<Haplotype>& haplotypes);
    void unprime() noexcept;
    bool is_primed() const noexcept;
    
    double evaluate(const Genotype<Haplotype>& genotype) const;
    double evaluate(const GenotypeIndex& genotype) const;
    
private:
    const HaplotypeLikelihoodCache& likelihoods_;
    std::vector<HaplotypeLikelihoodCache::LikelihoodVectorRef> indexed_likelihoods_;
    mutable std::vector<double> buffer_;
    
    // These are just for optimisation
    double evaluate_haploid(const Genotype<Haplotype>& genotype) const;
    double evaluate_diploid(const Genotype<Haplotype>& genotype) const;
    double evaluate_triploid(const Genotype<Haplotype>& genotype) const;
    double evaluate_tetraploid(const Genotype<Haplotype>& genotype) const;
    double evaluate_polyploid(const Genotype<Haplotype>& genotype) const;
};

template <typename G>
std::vector<double>&
evaluate(const std::vector<G>& genotypes, const GermlineLikelihoodModel& model, std::vector<double>& result)
{
    result.resize(genotypes.size());
    std::transform(std::cbegin(genotypes), std::cend(genotypes), std::begin(result),
                   [&] (const auto& genotype) { return model.evaluate(genotype); });
    return result;
}

template <typename G>
std::vector<double>
evaluate(const std::vector<G>& genotypes, const GermlineLikelihoodModel& model)
{
    std::vector<double> result {};
    evaluate(genotypes, model, result);
    return result;
}

} // namespace model
} // namespace octopus

#endif
