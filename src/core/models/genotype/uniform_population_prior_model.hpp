// Copyright (c) 2015-2018 Daniel Cooke
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
    using PopulationPriorModel::GenotypeReference;
    using PopulationPriorModel::GenotypeIndiceVectorReference;
    
    UniformPopulationPriorModel() = default;
    
    UniformPopulationPriorModel(const UniformPopulationPriorModel&)            = default;
    UniformPopulationPriorModel& operator=(const UniformPopulationPriorModel&) = default;
    UniformPopulationPriorModel(UniformPopulationPriorModel&&)                 = default;
    UniformPopulationPriorModel& operator=(UniformPopulationPriorModel&&)      = default;
    
    virtual ~UniformPopulationPriorModel() = default;

private:
    bool is_primed_ = false;
    
    LogProbability do_evaluate(const std::vector<Genotype<Haplotype>>& genotypes) const override
    {
        return 1.0;
    }
    LogProbability do_evaluate(const std::vector<GenotypeReference>& genotypes) const override
    {
        return 1.0;
    }
    LogProbability do_evaluate(const std::vector<GenotypeIndex>& genotypes) const override
    {
        return 1.0;
    }
    LogProbability do_evaluate(const std::vector<GenotypeIndiceVectorReference>& indices) const override
    {
        return 1.0;
    }
    void do_prime(const std::vector<Haplotype>& haplotypes) override
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
