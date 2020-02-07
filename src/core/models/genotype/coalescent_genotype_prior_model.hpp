// Copyright (c) 2015-2019 Daniel Cooke
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
    using GenotypePriorModel::LogProbability;
    
    CoalescentGenotypePriorModel() = delete;
    
    CoalescentGenotypePriorModel(CoalescentModel model) : model_ {std::move(model)} {}
    
    CoalescentGenotypePriorModel(const CoalescentGenotypePriorModel&)            = default;
    CoalescentGenotypePriorModel& operator=(const CoalescentGenotypePriorModel&) = default;
    CoalescentGenotypePriorModel(CoalescentGenotypePriorModel&&)                 = default;
    CoalescentGenotypePriorModel& operator=(CoalescentGenotypePriorModel&&)      = default;
    
    virtual ~CoalescentGenotypePriorModel() = default;

private:
    CoalescentModel model_;
    
    virtual LogProbability do_evaluate(const Genotype<Haplotype>& genotype) const override
    {
        return model_.evaluate(genotype);
    }
    virtual LogProbability do_evaluate(const Genotype<IndexedHaplotype<>>& genotype) const override
    {
        return model_.evaluate(genotype);
    }
    void do_prime(const HaplotypeBlock& haplotypes) override
    {
        model_.prime(haplotypes);
    }
    void do_unprime() noexcept override
    {
        model_.unprime();
    }
    bool check_is_primed() const noexcept override
    {
        return model_.is_primed();
    }
};

} // namespace octopus

#endif
