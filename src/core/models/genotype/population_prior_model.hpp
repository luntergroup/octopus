// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef population_prior_model_hpp
#define population_prior_model_hpp

#include <vector>
#include <functional>

#include "core/types/haplotype.hpp"

namespace octopus {

class PopulationPriorModel
{
public:
    PopulationPriorModel() = default;
    
    PopulationPriorModel(const PopulationPriorModel&)            = delete;
    PopulationPriorModel& operator=(const PopulationPriorModel&) = delete;
    PopulationPriorModel(PopulationPriorModel&&)                 = delete;
    PopulationPriorModel& operator=(PopulationPriorModel&&)      = delete;
    
    virtual ~PopulationPriorModel() = default;
    
    double evaluate(const std::vector<Haplotype>& haplotypes) const { return do_evaluate(haplotypes); }
    double evaluate(const std::vector<std::reference_wrapper<const Haplotype>>& haplotypes) const { return do_evaluate(haplotypes); }

private:
    virtual double do_evaluate(const std::vector<Haplotype>& haplotypes) const = 0;
    virtual double do_evaluate(const std::vector<std::reference_wrapper<const Haplotype>>& haplotypes) const = 0;
};

} // namespace octopus

#endif
