// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef hardy_weinberg_model_hpp
#define hardy_weinberg_model_hpp

#include <vector>
#include <unordered_map>
#include <functional>

#include <boost/optional.hpp>

#include "core/types/haplotype.hpp"
#include "core/types/indexed_haplotype.hpp"
#include "core/types/genotype.hpp"

namespace octopus {

class HardyWeinbergModel
{
public:
    using LogProbability = double;
    using GenotypeReference = std::reference_wrapper<const Genotype<IndexedHaplotype<>>>;
    using HaplotypeFrequencyVector = std::vector<double>;
    
    HardyWeinbergModel() = default;
    
    HardyWeinbergModel(IndexedHaplotype<> reference);
    HardyWeinbergModel(HaplotypeFrequencyVector haplotype_frequencies);
    
    HardyWeinbergModel(const HardyWeinbergModel&)            = default;
    HardyWeinbergModel& operator=(const HardyWeinbergModel&) = default;
    HardyWeinbergModel(HardyWeinbergModel&&)                 = default;
    HardyWeinbergModel& operator=(HardyWeinbergModel&&)      = default;
    
    ~HardyWeinbergModel() = default;
    
    void set_frequencies(HaplotypeFrequencyVector haplotype_frequencies);
    
    HaplotypeFrequencyVector& frequencies() noexcept;
    const HaplotypeFrequencyVector& frequencies() const noexcept;
    
    LogProbability evaluate(const Genotype<IndexedHaplotype<>>& genotype) const;
    LogProbability evaluate(const std::vector<Genotype<IndexedHaplotype<>>>& genotypes) const;
    LogProbability evaluate(const std::vector<GenotypeReference>& genotypes) const;
    
private:
    boost::optional<IndexedHaplotype<>> reference_;
    mutable HaplotypeFrequencyVector haplotype_frequencies_;
    mutable bool empirical_;
};

} // namespace octopus


#endif
