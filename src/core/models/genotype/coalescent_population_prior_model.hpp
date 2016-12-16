// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef coalescent_population_prior_model_hpp
#define coalescent_population_prior_model_hpp

#include <vector>
#include <utility>
#include <iterator>
#include <algorithm>
#include <functional>

#include "population_prior_model.hpp"
#include "../mutation/coalescent_model.hpp"

namespace octopus {

class CoalescentPopulationPriorModel : public PopulationPriorModel
{
public:
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
        thread_local std::vector<std::reference_wrapper<const Haplotype>> haplotypes {};
        haplotypes.clear();
        for (const auto& genotype : genotypes) {
            std::copy(std::cbegin(genotype), std::cend(genotype), std::back_inserter(haplotypes));
        }
        return model_.evaluate(haplotypes);
    }
    virtual double do_evaluate(const std::vector<std::reference_wrapper<const Genotype<Haplotype>>>& genotypes) const override
    {
        thread_local std::vector<std::reference_wrapper<const Haplotype>> haplotypes {};
        haplotypes.clear();
        for (const auto& genotype : genotypes) {
            std::copy(std::cbegin(genotype.get()), std::cend(genotype.get()), std::back_inserter(haplotypes));
        }
        return model_.evaluate(haplotypes);
    }
};
    
} // namespace octopus

#endif
