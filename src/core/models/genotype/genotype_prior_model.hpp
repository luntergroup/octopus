// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef genotype_prior_model_hpp
#define genotype_prior_model_hpp

#include "core/types/haplotype.hpp"
#include "core/types/genotype.hpp"

namespace octopus {

class GenotypePriorModel
{
public:
    GenotypePriorModel() = default;
    
    GenotypePriorModel(const GenotypePriorModel&)            = default;
    GenotypePriorModel& operator=(const GenotypePriorModel&) = default;
    GenotypePriorModel(GenotypePriorModel&&)                 = default;
    GenotypePriorModel& operator=(GenotypePriorModel&&)      = default;
    
    virtual ~GenotypePriorModel() = default;
    
    double evaluate(const Genotype<Haplotype>& genotype) const { return do_evaluate(genotype); }
    
private:
    virtual double do_evaluate(const Genotype<Haplotype>& genotype) const = 0;
};

} // namespace octopus

#endif
