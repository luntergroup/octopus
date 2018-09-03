// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef coalescent_population_prior_model_hpp
#define coalescent_population_prior_model_hpp

#include <vector>
#include <array>
#include <utility>
#include <iterator>
#include <algorithm>
#include <functional>

#include "population_prior_model.hpp"
#include "hardy_weinberg_model.hpp"
#include "../mutation/coalescent_model.hpp"

namespace octopus {

class CoalescentPopulationPriorModel : public PopulationPriorModel
{
public:
    using PopulationPriorModel::LogProbability;
    using PopulationPriorModel::GenotypeReference;
    using PopulationPriorModel::GenotypeIndiceVectorReference;
    
    CoalescentPopulationPriorModel() = delete;
    
    CoalescentPopulationPriorModel(CoalescentModel segregation_model);
    CoalescentPopulationPriorModel(CoalescentModel segregation_model, HardyWeinbergModel genotype_model);
    
    CoalescentPopulationPriorModel(const CoalescentPopulationPriorModel&)            = default;
    CoalescentPopulationPriorModel& operator=(const CoalescentPopulationPriorModel&) = default;
    CoalescentPopulationPriorModel(CoalescentPopulationPriorModel&&)                 = default;
    CoalescentPopulationPriorModel& operator=(CoalescentPopulationPriorModel&&)      = default;
    
    virtual ~CoalescentPopulationPriorModel() = default;

private:
    using HaplotypeReference = std::reference_wrapper<const Haplotype>;
    
    CoalescentModel segregation_model_;
    HardyWeinbergModel genotype_model_;
    
    mutable std::vector<unsigned> index_buffer_;
    
    LogProbability do_evaluate(const std::vector<Genotype<Haplotype>>& genotypes) const override
    {
        return evaluate_helper(genotypes);
    }
    LogProbability do_evaluate(const std::vector<GenotypeReference>& genotypes) const override
    {
        return evaluate_helper(genotypes);
    }
    LogProbability do_evaluate(const std::vector<GenotypeIndex>& indices) const override
    {
        return evaluate_helper(indices);
    }
    LogProbability do_evaluate(const std::vector<GenotypeIndiceVectorReference>& indices) const override
    {
        return evaluate_helper(indices);
    }
    void do_prime(const std::vector<Haplotype>& haplotypes) override
    {
        segregation_model_.prime(haplotypes);
    }
    void do_unprime() noexcept override
    {
        segregation_model_.unprime();
    }
    bool check_is_primed() const noexcept override
    {
        return segregation_model_.is_primed();
    }
    
    template <typename Range>
    LogProbability evaluate_helper(const Range& genotypes) const;
    template <typename Range>
    LogProbability evaluate_segregation_model(const Range& genotypes) const;
    LogProbability evaluate_segregation_model(const std::vector<std::vector<unsigned>>& indices) const;
    LogProbability evaluate_segregation_model(const std::vector<GenotypeIndiceVectorReference>& indices) const;
};

template <typename Range>
CoalescentPopulationPriorModel::LogProbability CoalescentPopulationPriorModel::evaluate_helper(const Range& genotypes) const
{
    // p({g_1, ..., g_n}) = p(g_1 u ... u g_n) p({g_1, ..., g_n} | g_1 u ... u g_n)
    // => ln p({g_1, ..., g_n}) = ln p(g_1 u ... u g_n) + ln p({g_1, ..., g_n} | g_1 u ... u g_n)
    // i.e The prior probability of observing a particular combination of genotypes is the
    // probability the haplotypes defined by the set of genotypes segregate, times the probability
    // of the particular genotypes given the haplotypes segregate.
    return evaluate_segregation_model(genotypes) + genotype_model_.evaluate(genotypes);
}

namespace detail {

template <typename Container>
void append(const Genotype<Haplotype>& genotype, Container& haplotypes)
{
    std::copy(std::cbegin(genotype), std::cend(genotype), std::back_inserter(haplotypes));
}

inline const Genotype<Haplotype>& get(const Genotype<Haplotype>& genotype) noexcept
{
    return genotype;
}

inline const Genotype<Haplotype>& get(const CoalescentPopulationPriorModel::GenotypeReference& genotype) noexcept
{
    return genotype.get();
}

inline const Haplotype& get(const Genotype<Haplotype>& genotype, const unsigned i) noexcept
{
    return genotype[i];
}

inline const Haplotype& get(const CoalescentPopulationPriorModel::GenotypeReference& genotype, const unsigned i) noexcept
{
    return genotype.get()[i];
}

inline auto ploidy(const Genotype<Haplotype>& genotype) noexcept
{
    return genotype.ploidy();
}

} // namespace detail

template <typename Range>
CoalescentPopulationPriorModel::LogProbability
CoalescentPopulationPriorModel::evaluate_segregation_model(const Range& genotypes) const
{
    if (genotypes.size() == 1) return segregation_model_.evaluate(detail::get(genotypes.front()));
    if (genotypes.size() == 2) {
        const auto ploidy1 = detail::ploidy(genotypes[0]);
        const auto ploidy2 = detail::ploidy(genotypes[1]);
        if (ploidy1 == ploidy2) {
            if (ploidy1 == 1) {
                using detail::get;
                const std::array<HaplotypeReference, 2> haplotypes {get(genotypes[0], 0), get(genotypes[0], 0)};
                return  segregation_model_.evaluate(haplotypes);
            } else if (ploidy1 == 2) {
                using detail::get;
                const std::array<HaplotypeReference, 4> haplotypes {get(genotypes[0], 0), get(genotypes[0], 1),
                                                                    get(genotypes[1], 0), get(genotypes[1], 1)};
                return segregation_model_.evaluate(haplotypes);
            }
        }
    }
    std::vector<HaplotypeReference> haplotypes {};
    haplotypes.reserve(10 * genotypes.size());
    for (const auto& genotype : genotypes) {
        detail::append(genotype, haplotypes);
    }
    return segregation_model_.evaluate(haplotypes);
}

} // namespace octopus

#endif
