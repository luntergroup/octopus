// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef population_prior_model_hpp
#define population_prior_model_hpp

#include <vector>
#include <functional>

#include "core/types/haplotype.hpp"
#include "core/types/genotype.hpp"

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
    
    double evaluate(const std::vector<Genotype<Haplotype>>& genotypes) const { return do_evaluate(genotypes); }
    double evaluate(const std::vector<std::reference_wrapper<const Genotype<Haplotype>>>& genotypes) const { return do_evaluate(genotypes); }

private:
    virtual double do_evaluate(const std::vector<Genotype<Haplotype>>& genotypes) const = 0;
    virtual double do_evaluate(const std::vector<std::reference_wrapper<const Genotype<Haplotype>>>& genotypes) const = 0;
};

} // namespace octopus

#endif
