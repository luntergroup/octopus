// Copyright (c) 2015-2021 Daniel Cooke
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
    using PopulationPriorModel::LogProbability;
    using PopulationPriorModel::GenotypeVector;
    using PopulationPriorModel::GenotypeReference;
    using PopulationPriorModel::GenotypeReferenceVector;
    
    UniformPopulationPriorModel() = default;
    
    UniformPopulationPriorModel(const UniformPopulationPriorModel&)            = delete;
    UniformPopulationPriorModel& operator=(const UniformPopulationPriorModel&) = delete;
    UniformPopulationPriorModel(UniformPopulationPriorModel&&)                 = delete;
    UniformPopulationPriorModel& operator=(UniformPopulationPriorModel&&)      = delete;
    
    virtual ~UniformPopulationPriorModel() = default;

private:
    bool is_primed_ = false;
    
    LogProbability do_evaluate(const GenotypeVector& genotypes) const override
    {
        return 1;
    }
    LogProbability do_evaluate(const GenotypeReferenceVector& genotypes) const override
    {
        return 1;
    }
    void do_prime(const HaplotypeBlock& haplotypes) override
    {
        is_primed_ = true;
    }
    void do_unprime() noexcept override
    {
        is_primed_ = false;
    }
    bool check_is_primed() const noexcept override
    {
        return is_primed_;
    }
};
    
} // namespace octopus

#endif
