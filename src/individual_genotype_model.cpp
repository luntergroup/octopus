//
//  individual_genotype_model.cpp
//  Octopus
//
//  Created by Daniel Cooke on 01/04/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#include "individual_genotype_model.hpp"

#include <utility>
#include <algorithm>
#include <cmath>
#include <cassert>
#include <iostream>

#include "fixed_ploidy_genotype_likelihood_model.hpp"
#include "maths.hpp"
#include "logging.hpp"

namespace Octopus
{
namespace GenotypeModel
{
    Individual::Individual(unsigned ploidy, const CoalescentModel& genotype_prior_model)
    :
    ploidy_ {ploidy},
    genotype_prior_model_ {genotype_prior_model}
    {}
    
    Individual::Latents::Latents(GenotypeProbabilityVector&& genotype_probabilities)
    :
    genotype_probabilities {std::move(genotype_probabilities)}
    {}
    
    Individual::InferredLatents::InferredLatents(Latents&& posteriors, double log_evidence)
    :
    posteriors {std::move(posteriors)},
    log_evidence {log_evidence}
    {}
    
    Individual::InferredLatents
    Individual::infer_latents(const SampleIdType& sample,
                              const std::vector<Genotype<Haplotype>>& genotypes,
                              const HaplotypeLikelihoodCache& haplotype_likelihoods) const
    {
        assert(!genotypes.empty());
        assert(genotypes.front().ploidy() == ploidy_);
        
        FixedPloidyGenotypeLikelihoodModel likelihood_model {ploidy_, haplotype_likelihoods};
        
        std::vector<double> result(genotypes.size());
        
        std::transform(std::cbegin(genotypes), std::cend(genotypes), std::begin(result),
                       [this, &sample, &likelihood_model] (const auto& genotype) {
                           return genotype_prior_model_.get().evaluate(genotype)
                                        + likelihood_model.log_likelihood(sample, genotype);
                       });
        
        auto log_evidence = Maths::normalise_exp(result);
        
        return InferredLatents {Latents {std::move(result)}, log_evidence};
    }
} // namesapce GenotypeModel
} // namespace Octopus
