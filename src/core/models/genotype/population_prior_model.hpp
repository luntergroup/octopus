// Copyright (c) 2015-2020 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef population_prior_model_hpp
#define population_prior_model_hpp

#include <vector>
#include <functional>
#include <initializer_list>

#include "core/types/haplotype.hpp"
#include "core/types/indexed_haplotype.hpp"
#include "core/types/genotype.hpp"
#include "containers/mappable_block.hpp"

namespace octopus {

class PopulationPriorModel
{
public:
    using LogProbability = double;
    using HaplotypeBlock = MappableBlock<Haplotype>;
    using GenotypeVector = std::vector<Genotype<IndexedHaplotype<>>>;
    using GenotypeReference = std::reference_wrapper<const Genotype<IndexedHaplotype<>>>;
    using GenotypeReferenceVector = std::vector<GenotypeReference>;
    
    PopulationPriorModel() = default;
    
    PopulationPriorModel(const PopulationPriorModel&)            = delete;
    PopulationPriorModel& operator=(const PopulationPriorModel&) = delete;
    PopulationPriorModel(PopulationPriorModel&&)                 = delete;
    PopulationPriorModel& operator=(PopulationPriorModel&&)      = delete;
    
    virtual ~PopulationPriorModel() = default;
    
    void prime(const HaplotypeBlock& haplotypes) { do_prime(haplotypes); }
    void unprime() noexcept { do_unprime(); }
    bool is_primed() const noexcept { return check_is_primed(); }
    
    LogProbability evaluate(const GenotypeVector& genotypes) const { return do_evaluate(genotypes); }
    LogProbability evaluate(const GenotypeReferenceVector& genotypes) const { return do_evaluate(genotypes); }
    
    LogProbability evaluate(std::initializer_list<GenotypeReference> genotypes) const
    {
        buffer_.assign(genotypes);
        return evaluate(buffer_);
    }
    
private:
    std::vector<Haplotype> haplotypes_;
    mutable std::vector<GenotypeReference> buffer_;
    
    virtual LogProbability do_evaluate(const GenotypeVector& genotypes) const
    {
        buffer_.clear();
        buffer_.reserve(genotypes.size());
        for (const auto& genotype : genotypes) {
            buffer_.push_back(std::cref(genotype));
        }
        return evaluate(buffer_);
    }
    virtual LogProbability do_evaluate(const GenotypeReferenceVector& genotypes) const = 0;
    virtual void do_prime(const HaplotypeBlock& haplotypes) {};
    virtual void do_unprime() noexcept {};
    virtual bool check_is_primed() const noexcept = 0;
};

} // namespace octopus

#endif
