// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef coalescent_population_prior_model_hpp
#define coalescent_population_prior_model_hpp

#include <vector>
#include <utility>
#include <iterator>
#include <algorithm>
#include <functional>

#include <boost/container/small_vector.hpp>

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

} // namespace detail

template <typename Container>
double CoalescentPopulationPriorModel::do_evaluate_helper(const Container& genotypes) const
{
    using HaplotypeReference = std::reference_wrapper<const Haplotype>;
    std::vector<HaplotypeReference> haplotypes {};
    haplotypes.reserve(genotypes.size() * 10);
    for (const auto& genotype : genotypes) {
        detail::append(genotype, haplotypes);
    }
    return model_.evaluate(haplotypes);
}

} // namespace octopus

#endif
