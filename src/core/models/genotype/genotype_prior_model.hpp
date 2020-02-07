// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef genotype_prior_model_hpp
#define genotype_prior_model_hpp

#include "core/types/haplotype.hpp"
#include "core/types/indexed_haplotype.hpp"
#include "core/types/genotype.hpp"
#include "containers/mappable_block.hpp"
#include "utils/maths.hpp"

namespace octopus {

class GenotypePriorModel
{
public:
    using LogProbability = double;
    using HaplotypeBlock = MappableBlock<Haplotype>;
    
    GenotypePriorModel() = default;
    
    GenotypePriorModel(const GenotypePriorModel&)            = default;
    GenotypePriorModel& operator=(const GenotypePriorModel&) = default;
    GenotypePriorModel(GenotypePriorModel&&)                 = default;
    GenotypePriorModel& operator=(GenotypePriorModel&&)      = default;
    
    virtual ~GenotypePriorModel() = default;
    
    void prime(const HaplotypeBlock& haplotypes) { do_prime(haplotypes); }
    void unprime() noexcept { do_unprime(); }
    bool is_primed() const noexcept { return check_is_primed(); }
    
    LogProbability evaluate(const Genotype<Haplotype>& genotype) const { return do_evaluate(genotype); }
    LogProbability evaluate(const Genotype<IndexedHaplotype<>>& genotype) const { return do_evaluate(genotype); }
    
private:
    virtual LogProbability do_evaluate(const Genotype<Haplotype>& genotype) const = 0;
    virtual LogProbability do_evaluate(const Genotype<IndexedHaplotype<>>& genotype) const = 0;
    virtual void do_prime(const HaplotypeBlock& haplotypes) {};
    virtual void do_unprime() noexcept {};
    virtual bool check_is_primed() const noexcept = 0;
};

// non-member methods

template <typename Container1, typename Container2>
Container2&
evaluate(const Container1& genotypes, const GenotypePriorModel& model, Container2& result,
         const bool normalise = false, const bool add = false)
{
    if (add) {
        assert(result.size() == genotypes.size());
        std::transform(std::cbegin(genotypes), std::cend(genotypes), std::cbegin(result), std::begin(result),
                       [&] (const auto& genotype, auto curr) { return curr + model.evaluate(genotype); });
    } else {
        result.resize(genotypes.size());
        std::transform(std::cbegin(genotypes), std::cend(genotypes), std::begin(result),
                       [&] (const auto& genotype) { return model.evaluate(genotype); });
    }
    if (normalise) maths::normalise_logs(result);
    return result;
}

template <typename Container>
auto evaluate(const Container& genotypes, const GenotypePriorModel& model, const bool normalise = false)
{
    std::vector<GenotypePriorModel::LogProbability> result(genotypes.size());
    evaluate(genotypes, model, result, normalise, false);
    return result;
}

} // namespace octopus

#endif
