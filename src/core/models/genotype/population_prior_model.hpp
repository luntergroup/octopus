// Copyright (c) 2016 Daniel Cooke
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
    
    double evaluate(const std::vector<Genotype<Haplotype>>& genotypes) const { return do_evaluate(genotypes); }
    double evaluate(const std::vector<GenotypeReference>& genotypes) const { return do_evaluate(genotypes); }
    double evaluate(const std::vector<std::vector<unsigned>>& indices) const { return do_evaluate(indices); }
    double evaluate(const std::vector<GenotypeIndiceVectorReference>& indices) const { return do_evaluate(indices); }
    
private:
    std::vector<Haplotype> haplotypes_;
    
    virtual double do_evaluate(const std::vector<Genotype<Haplotype>>& genotypes) const = 0;
    virtual double do_evaluate(const std::vector<GenotypeReference>& genotypes) const = 0;
    virtual double do_evaluate(const std::vector<std::vector<unsigned>>& indices) const = 0;
    virtual double do_evaluate(const std::vector<GenotypeIndiceVectorReference>& indices) const = 0;
    virtual void do_prime(const std::vector<Haplotype>& haplotypes) {};
    virtual void do_unprime() noexcept {};
    virtual bool check_is_primed() const noexcept = 0;
};

} // namespace octopus

#endif
