// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef hardy_weinberg_model_hpp
#define hardy_weinberg_model_hpp

#include <vector>
#include <unordered_map>
#include <functional>

#include <boost/optional.hpp>

#include "core/types/haplotype.hpp"
#include "core/types/genotype.hpp"

namespace octopus {

class HardyWeinbergModel
{
public:
    using LogProbability               = double;
    using GenotypeReference            = std::reference_wrapper<const Genotype<Haplotype>>;
    using GenotypeReferenceVector      = std::vector<GenotypeReference>;
    using GenotypeIndexReference       = std::reference_wrapper<const GenotypeIndex>;
    using GenotypeIndexVector          = std::vector<GenotypeIndex>;
    using GenotypeIndexReferenceVector = std::vector<GenotypeIndexReference>;
    
    using HaplotypeFrequencyMap    = std::unordered_map<Haplotype, double>;
    using HaplotypeFrequencyVector = std::vector<double>;
    
    HardyWeinbergModel() = default;
    
    HardyWeinbergModel(Haplotype reference);
    HardyWeinbergModel(unsigned reference_idx);
    HardyWeinbergModel(HaplotypeFrequencyMap haplotype_frequencies);
    HardyWeinbergModel(HaplotypeFrequencyVector haplotype_frequencies);
    
    HardyWeinbergModel(const HardyWeinbergModel&)            = default;
    HardyWeinbergModel& operator=(const HardyWeinbergModel&) = default;
    HardyWeinbergModel(HardyWeinbergModel&&)                 = default;
    HardyWeinbergModel& operator=(HardyWeinbergModel&&)      = default;
    
    ~HardyWeinbergModel() = default;
    
    void set_frequencies(HaplotypeFrequencyMap haplotype_frequencies);
    void set_frequencies(HaplotypeFrequencyVector haplotype_frequencies);
    
    HaplotypeFrequencyMap& frequencies() noexcept;
    HaplotypeFrequencyVector& index_frequencies() noexcept;
    
    LogProbability evaluate(const Genotype<Haplotype>& genotype) const;
    LogProbability evaluate(const GenotypeIndex& genotype) const;
    
    LogProbability evaluate(const std::vector<Genotype<Haplotype>>& genotypes) const;
    LogProbability evaluate(const GenotypeReferenceVector& genotypes) const;
    LogProbability evaluate(const GenotypeIndexVector& genotypes) const;
    LogProbability evaluate(const GenotypeIndexReferenceVector& genotypes) const;
    
private:
    boost::optional<Haplotype> reference_;
    boost::optional<unsigned> reference_idx_;
    mutable HaplotypeFrequencyMap haplotype_frequencies_;
    mutable HaplotypeFrequencyVector haplotype_idx_frequencies_;
    mutable bool empirical_;
};

} // namespace octopus


#endif
