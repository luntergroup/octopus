// Copyright (c) 2015-2020 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "coalescent_population_prior_model.hpp"

namespace octopus {

CoalescentPopulationPriorModel::CoalescentPopulationPriorModel(CoalescentModel segregation_model)
: segregation_model_ {std::move(segregation_model)}
, genotype_model_ {}
{}

CoalescentPopulationPriorModel::CoalescentPopulationPriorModel(CoalescentModel segregation_model, HardyWeinbergModel genotype_model)
: segregation_model_ {std::move(segregation_model)}
, genotype_model_ {std::move(genotype_model)}
{}

CoalescentPopulationPriorModel::LogProbability
CoalescentPopulationPriorModel::do_evaluate(const GenotypeReferenceVector& genotypes) const
{
    // p({g_1, ..., g_n}) = p(g_1 u ... u g_n) p({g_1, ..., g_n} | g_1 u ... u g_n)
    // => ln p({g_1, ..., g_n}) = ln p(g_1 u ... u g_n) + ln p({g_1, ..., g_n} | g_1 u ... u g_n)
    // i.e The prior probability of observing a particular combination of genotypes is the
    // probability the haplotypes defined by the set of genotypes segregate, times the probability
    // of the particular genotypes given the haplotypes segregate.
    return evaluate_segregation_model(genotypes) + genotype_model_.evaluate(genotypes);
}

template <typename Range>
auto sum_ploidies(const Range& genotypes) noexcept
{
    return std::accumulate(std::cbegin(genotypes), std::cend(genotypes), 0u,
                           [] (auto total, const auto& genotype) noexcept { return total + ploidy(genotype.get()); });
}

CoalescentPopulationPriorModel::LogProbability
CoalescentPopulationPriorModel::evaluate_segregation_model(const GenotypeReferenceVector& genotypes) const
{
    if (genotypes.size() == 1) {
        return segregation_model_.evaluate(genotypes.front().get());
    }
    const auto total_ploidy = sum_ploidies(genotypes);
    haplotype_buffer_.clear();
    haplotype_buffer_.reserve(total_ploidy);
    for (const auto& genotype : genotypes) {
        for (const auto& haplotype : genotype.get()) {
            haplotype_buffer_.push_back(haplotype);
        }
    }
    return segregation_model_.evaluate(haplotype_buffer_);
}

} // namespace
