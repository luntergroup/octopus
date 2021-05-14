// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef pedigree_prior_model_hpp
#define pedigree_prior_model_hpp

#include <vector>
#include <memory>

#include "basics/pedigree.hpp"
#include "population_prior_model.hpp"
#include "inheritance_model.hpp"

namespace octopus {

class PedigreePriorModel
{
public:
    using SampleVector   = std::vector<SampleName>;
    using LogProbability = double;

    PedigreePriorModel() = delete;

    PedigreePriorModel(const SampleVector& samples,
                       const Pedigree& pedigree,
                       std::unique_ptr<PopulationPriorModel> population_prior_model,
                       InheritanceModel inheritance_model);
    
    PedigreePriorModel(const PedigreePriorModel&)            = delete;
    PedigreePriorModel& operator=(const PedigreePriorModel&) = delete;
    PedigreePriorModel(PedigreePriorModel&&)                 = default;
    PedigreePriorModel& operator=(PedigreePriorModel&&)      = default;
    
    ~PedigreePriorModel() = default;
    
    const PopulationPriorModel& population_model() const noexcept;
    const InheritanceModel& inheritance_model() const noexcept;

    void prime(const MappableBlock<Haplotype>& haplotypes);
    void unprime() noexcept;
    bool is_primed() const noexcept;

    LogProbability evaluate(const std::vector<Genotype<IndexedHaplotype<>>>& genotypes) const;

private:
    std::vector<std::size_t> founders_;
    std::vector<std::pair<std::size_t, std::pair<std::size_t, std::size_t>>> offspring_with_two_parents_;
    std::vector<std::pair<std::size_t, std::size_t>> offspring_with_one_parent_;
    std::vector<std::pair<std::pair<std::size_t, std::size_t>, std::pair<std::size_t, std::size_t>>> twins_with_two_parents_;
    std::vector<std::pair<std::pair<std::size_t, std::size_t>, std::size_t>> twins_with_one_parent_;
    std::unique_ptr<PopulationPriorModel> population_prior_model_;
    InheritanceModel inheritance_model_;
};

} // namespace octopus

#endif
