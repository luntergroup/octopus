// Copyright (c) 2016 Daniel Cooke
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
#include "../mutation/coalescent_model.hpp"

#include "timers.hpp"

namespace octopus {

class CoalescentPopulationPriorModel : public PopulationPriorModel
{
public:
    using PopulationPriorModel::GenotypeReference;
    
    CoalescentPopulationPriorModel() = delete;
    
    CoalescentPopulationPriorModel(CoalescentModel model) : model_ {std::move(model)} {}
    
    CoalescentPopulationPriorModel(const CoalescentPopulationPriorModel&)            = default;
    CoalescentPopulationPriorModel& operator=(const CoalescentPopulationPriorModel&) = default;
    CoalescentPopulationPriorModel(CoalescentPopulationPriorModel&&)                 = default;
    CoalescentPopulationPriorModel& operator=(CoalescentPopulationPriorModel&&)      = default;
    
    virtual ~CoalescentPopulationPriorModel() = default;

private:
    using HaplotypeReference = std::reference_wrapper<const Haplotype>;
    
    CoalescentModel model_;
        
    virtual double do_evaluate(const std::vector<Genotype<Haplotype>>& genotypes) const override
    {
        return do_evaluate_helper(genotypes);
    }
    virtual double do_evaluate(const std::vector<GenotypeReference>& genotypes) const override
    {
        return do_evaluate_helper(genotypes);
    }
    
    template <typename Container>
    double do_evaluate_helper(const Container& genotypes) const;
    
};

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

template <typename Container>
double CoalescentPopulationPriorModel::do_evaluate_helper(const Container& genotypes) const
{
    if (genotypes.size() == 1) return model_.evaluate(detail::get(genotypes.front()));
    if (genotypes.size() == 2) {
        const auto ploidy1 = detail::ploidy(genotypes[0]);
        const auto ploidy2 = detail::ploidy(genotypes[1]);
        if (ploidy1 == ploidy2) {
            if (ploidy1 == 1) {
                using detail::get;
                const std::array<HaplotypeReference, 2> haplotypes {get(genotypes[0], 0), get(genotypes[0], 0)};
                return  model_.evaluate(haplotypes);
            } else if (ploidy1 == 2) {
                using detail::get;
                const std::array<HaplotypeReference, 4> haplotypes {get(genotypes[0], 0), get(genotypes[0], 1),
                                                                    get(genotypes[1], 0), get(genotypes[1], 1)};
                return model_.evaluate(haplotypes);
            }
        }
    }
    std::vector<HaplotypeReference> haplotypes {};
    haplotypes.reserve(10 * genotypes.size());
    for (const auto& genotype : genotypes) {
        detail::append(genotype, haplotypes);
    }
    return model_.evaluate(haplotypes);
}

} // namespace octopus

#endif
