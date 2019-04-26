// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef somatic_mutation_model_hpp
#define somatic_mutation_model_hpp

#include <type_traits>
#include <functional>

#include "core/types/haplotype.hpp"
#include "containers/mappable_block.hpp"
#include "denovo_model.hpp"

namespace octopus {

class SomaticMutationModel
{
public:
    using LogProbability  = DeNovoModel::LogProbability;
    using Parameters      = DeNovoModel::Parameters;
    using CachingStrategy = DeNovoModel::CachingStrategy;
    
    SomaticMutationModel() = delete;
    
    SomaticMutationModel(Parameters params,
                         std::size_t num_haplotypes_hint = 1000,
                         CachingStrategy caching = CachingStrategy::value);
    
    SomaticMutationModel(const SomaticMutationModel&)            = default;
    SomaticMutationModel& operator=(const SomaticMutationModel&) = default;
    SomaticMutationModel(SomaticMutationModel&&)                 = default;
    SomaticMutationModel& operator=(SomaticMutationModel&&)      = default;
    
    ~SomaticMutationModel() = default;
    
    void prime(MappableBlock<Haplotype> haplotypes);
    void unprime() noexcept;
    bool is_primed() const noexcept;
    
    // ln p(somatic | germline)
    LogProbability evaluate(const Haplotype& somatic, const Haplotype& germline) const;
    LogProbability evaluate(unsigned somatic, unsigned germline) const;
    
private:
    DeNovoModel model_;
};

} // namespace octopus

#endif
