// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef inheritance_model_hpp
#define inheritance_model_hpp

#include "core/types/haplotype.hpp"
#include "core/types/indexed_haplotype.hpp"
#include "core/types/genotype.hpp"
#include "core/models/mutation/denovo_model.hpp"

namespace octopus {

class InheritanceModel
{
public:
    using LogProbability = double;
    
    InheritanceModel() = delete;
    
    InheritanceModel(DeNovoModel mutation_model);
    
    InheritanceModel(const InheritanceModel&)            = default;
    InheritanceModel& operator=(const InheritanceModel&) = default;
    InheritanceModel(InheritanceModel&&)                 = default;
    InheritanceModel& operator=(InheritanceModel&&)      = default;
    
    ~InheritanceModel() = default;
    
    void prime(const MappableBlock<Haplotype>& haplotypes);
    void unprime() noexcept;
    bool is_primed() const noexcept;
    
    // p(offspring | parent)
    LogProbability evaluate(const Genotype<Haplotype>& offspring, const Genotype<Haplotype>& parent) const;
    LogProbability evaluate(const Genotype<IndexedHaplotype<>>& offspring, const Genotype<IndexedHaplotype<>>& parent) const;
    
    // p(offspring | mother, father)
    LogProbability evaluate(const Genotype<Haplotype>& offspring, const Genotype<Haplotype>& mother, const Genotype<Haplotype>& father) const;
    LogProbability evaluate(const Genotype<IndexedHaplotype<>>& offspring, const Genotype<IndexedHaplotype<>>& mother, const Genotype<IndexedHaplotype<>>& father) const;

    // p(twin1, twin2 | parent)
    LogProbability evaluate_twins(const Genotype<Haplotype>& twin1, const Genotype<Haplotype>& twin2, const Genotype<Haplotype>& parent) const;
    LogProbability evaluate_twins(const Genotype<IndexedHaplotype<>>& twin1, const Genotype<IndexedHaplotype<>>& twin2, const Genotype<IndexedHaplotype<>>& parent) const;
    
    // p(twin1, twin2 | mother, father)
    LogProbability evaluate_twins(const Genotype<Haplotype>& twin1, const Genotype<Haplotype>& twin2, const Genotype<Haplotype>& mother, const Genotype<Haplotype>& father) const;
    LogProbability evaluate_twins(const Genotype<IndexedHaplotype<>>& twin1, const Genotype<IndexedHaplotype<>>& twin2, const Genotype<IndexedHaplotype<>>& mother, const Genotype<IndexedHaplotype<>>& father) const;

private:
    DeNovoModel mutation_model_;
};

} // namespace octopus

#endif
