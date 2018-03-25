// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef genotype_prior_model_hpp
#define genotype_prior_model_hpp

#include "core/types/haplotype.hpp"
#include "core/types/genotype.hpp"

namespace octopus {

class GenotypePriorModel
{
public:
    GenotypePriorModel() = default;
    
    GenotypePriorModel(const GenotypePriorModel&)            = default;
    GenotypePriorModel& operator=(const GenotypePriorModel&) = default;
    GenotypePriorModel(GenotypePriorModel&&)                 = default;
    GenotypePriorModel& operator=(GenotypePriorModel&&)      = default;
    
    virtual ~GenotypePriorModel() = default;
    
    void prime(const std::vector<Haplotype>& haplotypes) { do_prime(haplotypes); }
    void unprime() noexcept { do_unprime(); }
    bool is_primed() const noexcept { return check_is_primed(); }
    
    double evaluate(const Genotype<Haplotype>& genotype) const { return do_evaluate(genotype); }
    double evaluate(const std::vector<unsigned>& genotype_indices) const { return do_evaluate(genotype_indices); }
    
private:
    virtual double do_evaluate(const Genotype<Haplotype>& genotype) const = 0;
    virtual double do_evaluate(const std::vector<unsigned>& genotype) const = 0;
    virtual void do_prime(const std::vector<Haplotype>& haplotypes) {};
    virtual void do_unprime() noexcept {};
    virtual bool check_is_primed() const noexcept = 0;
};

} // namespace octopus

#endif
