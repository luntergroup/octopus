//
//  known_ploidy_genotype_likelihood_model.hpp
//  Octopus
//
//  Created by Daniel Cooke on 01/04/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__known_ploidy_genotype_likelihood_model__
#define __Octopus__known_ploidy_genotype_likelihood_model__

#include <cstddef>
#include <iterator>
#include <vector>

#include "haplotype.hpp"
#include "genotype.hpp"
#include "aligned_read.hpp"
#include "haplotype_likelihood_cache.hpp"
#include "probability_matrix.hpp"

namespace octopus
{
namespace model
{
    class KnownPloidyGenotypeLikelihoodModel
    {
    public:
        KnownPloidyGenotypeLikelihoodModel() = delete;
        
        KnownPloidyGenotypeLikelihoodModel(unsigned ploidy, const HaplotypeLikelihoodCache& likelihoods);
        
        KnownPloidyGenotypeLikelihoodModel(const KnownPloidyGenotypeLikelihoodModel&)            = default;
        KnownPloidyGenotypeLikelihoodModel& operator=(const KnownPloidyGenotypeLikelihoodModel&) = default;
        KnownPloidyGenotypeLikelihoodModel(KnownPloidyGenotypeLikelihoodModel&&)                 = default;
        KnownPloidyGenotypeLikelihoodModel& operator=(KnownPloidyGenotypeLikelihoodModel&&)      = default;
        
        ~KnownPloidyGenotypeLikelihoodModel() = default;
        
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

#endif /* defined(__Octopus__known_ploidy_genotype_likelihood_model__) */
