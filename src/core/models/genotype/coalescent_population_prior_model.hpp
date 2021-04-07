// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef coalescent_population_prior_model_hpp
#define coalescent_population_prior_model_hpp

#include <vector>
#include <array>
#include <utility>
#include <iterator>
#include <algorithm>
#include <functional>

#include "population_prior_model.hpp"
#include "hardy_weinberg_model.hpp"
#include "../mutation/coalescent_model.hpp"

namespace octopus {

class CoalescentPopulationPriorModel : public PopulationPriorModel
{
public:
    using PopulationPriorModel::LogProbability;
    using PopulationPriorModel::GenotypeReference;
    using PopulationPriorModel::GenotypeReferenceVector;
    
    CoalescentPopulationPriorModel() = delete;
    
    CoalescentPopulationPriorModel(CoalescentModel segregation_model);
    CoalescentPopulationPriorModel(CoalescentModel segregation_model, HardyWeinbergModel genotype_model);
    
    CoalescentPopulationPriorModel(const CoalescentPopulationPriorModel&)            = delete;
    CoalescentPopulationPriorModel& operator=(const CoalescentPopulationPriorModel&) = delete;
    CoalescentPopulationPriorModel(CoalescentPopulationPriorModel&&)                 = delete;
    CoalescentPopulationPriorModel& operator=(CoalescentPopulationPriorModel&&)      = delete;
    
    virtual ~CoalescentPopulationPriorModel() = default;

private:
    CoalescentModel segregation_model_;
    HardyWeinbergModel genotype_model_;
    
    mutable std::vector<IndexedHaplotype<>> haplotype_buffer_;
    
    LogProbability do_evaluate(const GenotypeReferenceVector& genotypes) const override;
    
    void do_prime(const HaplotypeBlock& haplotypes) override
    {
        segregation_model_.prime(haplotypes);
    }
    void do_unprime() noexcept override
    {
        segregation_model_.unprime();
    }
    bool check_is_primed() const noexcept override
    {
        return segregation_model_.is_primed();
    }
    
    LogProbability evaluate_segregation_model(const GenotypeReferenceVector& genotypes) const;
};

} // namespace octopus

#endif
