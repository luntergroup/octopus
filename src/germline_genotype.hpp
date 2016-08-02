//
//  germline_genotype.hpp
//  Octopus
//
//  Created by Daniel Cooke on 01/04/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__germline_genotype__
#define __Octopus__germline_genotype__

#include "haplotype.hpp"
#include "genotype.hpp"
#include "haplotype_likelihood_cache.hpp"

namespace octopus { namespace model {

class GermlineGenotype
{
public:
    GermlineGenotype() = delete;
    
    GermlineGenotype(unsigned ploidy, const HaplotypeLikelihoodCache& likelihoods);
    
    GermlineGenotype(const GermlineGenotype&)            = default;
    GermlineGenotype& operator=(const GermlineGenotype&) = default;
    GermlineGenotype(GermlineGenotype&&)                 = default;
    GermlineGenotype& operator=(GermlineGenotype&&)      = default;
    
    ~GermlineGenotype() = default;
    
    double ln_likelihood(const Genotype<Haplotype>& genotype) const;
    
private:
    const HaplotypeLikelihoodCache& likelihoods_;
    
    unsigned ploidy_;
    
    // These are just for optimisation
    double ln_likelihood_haploid(const Genotype<Haplotype>& genotype) const;
    double ln_likelihood_diploid(const Genotype<Haplotype>& genotype) const;
    double ln_likelihood_triploid(const Genotype<Haplotype>& genotype) const;
    double ln_likelihood_tetraploid(const Genotype<Haplotype>& genotype) const;
    double ln_likelihood_polyploid(const Genotype<Haplotype>& genotype) const;
};

} // namespace model
} // namespace octopus

#endif /* defined(__Octopus__germline_genotype__) */
