// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef coalescent_genotype_prior_model_hpp
#define coalescent_genotype_prior_model_hpp

#include <utility>

#include "genotype_prior_model.hpp"
#include "../mutation/coalescent_model.hpp"

namespace octopus {

class CoalescentGenotypePriorModel : public GenotypePriorModel
{
public:
    CoalescentGenotypePriorModel() = delete;
    
    CoalescentGenotypePriorModel(CoalescentModel model) : model_ {std::move(model)} {}
    
    CoalescentGenotypePriorModel(const CoalescentGenotypePriorModel&)            = default;
    CoalescentGenotypePriorModel& operator=(const CoalescentGenotypePriorModel&) = default;
    CoalescentGenotypePriorModel(CoalescentGenotypePriorModel&&)                 = default;
    CoalescentGenotypePriorModel& operator=(CoalescentGenotypePriorModel&&)      = default;
    
    virtual ~CoalescentGenotypePriorModel() = default;

private:
    CoalescentModel model_;
    
    virtual double do_evaluate(const Genotype<Haplotype>& genotype) const override
    {
        return model_.evaluate(genotype);
    }
};

} // namespace octopus

#endif
