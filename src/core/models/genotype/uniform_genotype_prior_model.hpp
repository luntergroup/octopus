// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef uniform_genotype_prior_model_hpp
#define uniform_genotype_prior_model_hpp

#include "genotype_prior_model.hpp"

namespace octopus {

class UniformGenotypePriorModel : public GenotypePriorModel
{
public:
    UniformGenotypePriorModel() = default;
    
    UniformGenotypePriorModel(const UniformGenotypePriorModel&)            = default;
    UniformGenotypePriorModel& operator=(const UniformGenotypePriorModel&) = default;
    UniformGenotypePriorModel(UniformGenotypePriorModel&&)                 = default;
    UniformGenotypePriorModel& operator=(UniformGenotypePriorModel&&)      = default;
    
    virtual ~UniformGenotypePriorModel() = default;
    
private:
    virtual double do_evaluate(const Genotype<Haplotype>& genotype) const override { return 1.0; }
};

} // namespace octopus

#endif
