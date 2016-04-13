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
    
    Individual::Latents::Latents(const GenotypeVector& genotypes,
                                 GenotypeProbabilityVector&& genotype_posteriors)
    :
    genotypes {std::cref(genotypes)},
    genotype_posteriors {std::move(genotype_posteriors)}
    {}
    
    Individual::Latents
    Individual::infer_latents(const SampleIdType& sample,
                              const std::vector<Genotype<Haplotype>>& genotypes,
                              const HaplotypeLikelihoodCache& haplotype_likelihoods) const
    {
        assert(!genotypes.empty());
        
        FixedPloidyGenotypeLikelihoodModel likelihood_model {ploidy_, haplotype_likelihoods};
        
        std::vector<double> result(genotypes.size());
        
        std::transform(std::cbegin(genotypes), std::cend(genotypes), std::begin(result),
                       [this, &sample, &likelihood_model] (const auto& genotype) {
                           return std::log(genotype_prior_model_.get().evaluate(genotype))
                                        + likelihood_model.log_likelihood(sample, genotype);
                       });
        
        const auto log_evidence = Maths::log_sum_exp(result);
        
        //std::cout << "log_evidence = " << log_evidence << std::endl;
        
        Maths::normalise_exp(result);
        
        return Latents {genotypes, std::move(result)};
    }
    
    double Individual::log_evidence(const SampleIdType& sample,
                                    const HaplotypeLikelihoodCache& haplotype_likelihoods,
                                    const Latents& latents) const
    {
        FixedPloidyGenotypeLikelihoodModel likelihood_model {ploidy_, haplotype_likelihoods};
        
        std::vector<double> log_jps(latents.genotype_posteriors.size());
        
        std::transform(std::cbegin(latents.genotypes.get()), std::cend(latents.genotypes.get()),
                       std::cbegin(latents.genotype_posteriors),
                       std::begin(log_jps),
                       [&sample, &likelihood_model] (const auto& g, const auto& p) {
                           return std::log(p) + likelihood_model.log_likelihood(sample, g);
                       });
        
        return Maths::log_sum_exp(log_jps);
    }
} // namesapce GenotypeModel
} // namespace Octopus
