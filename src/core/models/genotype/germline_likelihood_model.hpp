// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef germline_likelihood_model_hpp
#define germline_likelihood_model_hpp

#include "core/types/haplotype.hpp"
#include "core/types/genotype.hpp"
#include "core/models/haplotype_likelihood_cache.hpp"

namespace octopus { namespace model {

class GermlineLikelihoodModel
{
public:
    GermlineLikelihoodModel() = delete;
    
    GermlineLikelihoodModel(const HaplotypeLikelihoodCache& likelihoods);
    
    GermlineLikelihoodModel(const GermlineLikelihoodModel&)            = default;
    GermlineLikelihoodModel& operator=(const GermlineLikelihoodModel&) = default;
    GermlineLikelihoodModel(GermlineLikelihoodModel&&)                 = default;
    GermlineLikelihoodModel& operator=(GermlineLikelihoodModel&&)      = default;
    
    ~GermlineLikelihoodModel() = default;
    
    double evaluate(const Genotype<Haplotype>& genotype) const;
    
private:
    const HaplotypeLikelihoodCache& likelihoods_;
    
    // These are just for optimisation
    double evaluate_haploid(const Genotype<Haplotype>& genotype) const;
    double evaluate_diploid(const Genotype<Haplotype>& genotype) const;
    double evaluate_triploid(const Genotype<Haplotype>& genotype) const;
    double evaluate_tetraploid(const Genotype<Haplotype>& genotype) const;
    double evaluate_polyploid(const Genotype<Haplotype>& genotype) const;
};

} // namespace model
} // namespace octopus

#endif
