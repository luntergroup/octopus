// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef uniform_population_prior_model_hpp
#define uniform_population_prior_model_hpp

#include <vector>
#include <utility>

#include "population_prior_model.hpp"
#include "../mutation/coalescent_model.hpp"

namespace octopus {

class UniformPopulationPriorModel : public PopulationPriorModel
{
public:
    UniformPopulationPriorModel() = default;
    
    UniformPopulationPriorModel(const UniformPopulationPriorModel&)            = default;
    UniformPopulationPriorModel& operator=(const UniformPopulationPriorModel&) = default;
    UniformPopulationPriorModel(UniformPopulationPriorModel&&)                 = default;
    UniformPopulationPriorModel& operator=(UniformPopulationPriorModel&&)      = default;
    
    virtual ~UniformPopulationPriorModel() = default;

private:
    virtual double do_evaluate(const std::vector<Haplotype>& haplotypes) const override
    {
        return 1.0;
    }
    virtual double do_evaluate(const std::vector<std::reference_wrapper<const Haplotype>>& haplotypes) const override
    {
        return 1.0;
    }
};
    
} // namespace octopus

#endif
