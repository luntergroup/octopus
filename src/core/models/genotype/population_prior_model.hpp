// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef population_prior_model_hpp
#define population_prior_model_hpp

#include <vector>
#include <functional>

#include "core/types/haplotype.hpp"
#include "core/types/genotype.hpp"

namespace octopus {

class PopulationPriorModel
{
public:
    using LogProbability = double;
    using GenotypeReference = std::reference_wrapper<const Genotype<Haplotype>>;
    using GenotypeIndiceVectorReference = std::reference_wrapper<const std::vector<unsigned>>;
    
    PopulationPriorModel() = default;
    
    PopulationPriorModel(const PopulationPriorModel&)            = delete;
    PopulationPriorModel& operator=(const PopulationPriorModel&) = delete;
    PopulationPriorModel(PopulationPriorModel&&)                 = delete;
    PopulationPriorModel& operator=(PopulationPriorModel&&)      = delete;
    
    virtual ~PopulationPriorModel() = default;
    
    void prime(const std::vector<Haplotype>& haplotypes) { do_prime(haplotypes); }
    void unprime() noexcept { do_unprime(); }
    bool is_primed() const noexcept { return check_is_primed(); }
    
    LogProbability evaluate(const std::vector<Genotype<Haplotype>>& genotypes) const { return do_evaluate(genotypes); }
    LogProbability evaluate(const std::vector<GenotypeReference>& genotypes) const { return do_evaluate(genotypes); }
    LogProbability evaluate(const std::vector<GenotypeIndex>& indices) const { return do_evaluate(indices); }
    LogProbability evaluate(const std::vector<GenotypeIndiceVectorReference>& genotypes) const { return do_evaluate(genotypes); }
    
private:
    std::vector<Haplotype> haplotypes_;
    
    virtual LogProbability do_evaluate(const std::vector<Genotype<Haplotype>>& genotypes) const = 0;
    virtual LogProbability do_evaluate(const std::vector<GenotypeReference>& genotypes) const = 0;
    virtual LogProbability do_evaluate(const std::vector<GenotypeIndex>& genotypes) const = 0;
    virtual LogProbability do_evaluate(const std::vector<GenotypeIndiceVectorReference>& indices) const = 0;
    virtual void do_prime(const std::vector<Haplotype>& haplotypes) {};
    virtual void do_unprime() noexcept {};
    virtual bool check_is_primed() const noexcept = 0;
};

} // namespace octopus

#endif
