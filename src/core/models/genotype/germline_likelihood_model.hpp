// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef __Octopus__germline_likelihood_model__
#define __Octopus__germline_likelihood_model__

#include <core/types/haplotype.hpp>
#include <core/types/genotype.hpp>
#include <core/models/haplotype_likelihood_cache.hpp>

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
    
    double ln_likelihood(const Genotype<Haplotype>& genotype) const;
    
private:
    const HaplotypeLikelihoodCache& likelihoods_;
    
    // These are just for optimisation
    double ln_likelihood_haploid(const Genotype<Haplotype>& genotype) const;
    double ln_likelihood_diploid(const Genotype<Haplotype>& genotype) const;
    double ln_likelihood_triploid(const Genotype<Haplotype>& genotype) const;
    double ln_likelihood_tetraploid(const Genotype<Haplotype>& genotype) const;
    double ln_likelihood_polyploid(const Genotype<Haplotype>& genotype) const;
};

} // namespace model
} // namespace octopus

#endif /* defined(__Octopus__germline_likelihood_model__) */
