//
//  individual_genotype_model.cpp
//  Octopus
//
//  Created by Daniel Cooke on 01/04/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#include "individual_genotype_model.hpp"

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
        Individual::Individual(unsigned ploidy)
        :
        ploidy_ {ploidy}
        {}
        
        auto calculate_haplotype_posteriors(const std::vector<Haplotype>& haplotypes,
                                            const std::vector<Genotype<Haplotype>>& genotypes,
                                            const std::vector<double>& genotype_log_posteriors)
        {
            std::unordered_map<std::reference_wrapper<const Haplotype>, double> result {haplotypes.size()};
            
            for (const auto& haplotype : haplotypes) {
                result.emplace(std::cref(haplotype), 0.0);
            }
            
            auto it = std::cbegin(genotype_log_posteriors);
            
            for (const auto& genotype : genotypes) {
                for (const auto& haplotype : genotype.copy_unique_ref()) {
                    result.at(haplotype) += *it;
                }
                ++it;
            }
            
            return result;
        }
        
        Individual::Latents
        Individual::infer_latents(const SampleIdType& sample,
                                  const std::vector<Haplotype>& haplotypes,
                                  const CoalescentModel& haplotype_model,
                                  const HaplotypeLikelihoodCache& haplotype_likelihoods) const
        {
            assert(!haplotypes.empty());
            
            auto genotypes = generate_all_genotypes(haplotypes, ploidy_);
            
            Logging::DebugLogger log {};
            if (DEBUG_MODE) {
                stream(log) << "There are " << genotypes.size() << " genotypes";
            }
            
            assert(!genotypes.empty());
            
            FixedPloidyGenotypeLikelihoodModel likelihood_model {ploidy_, haplotype_likelihoods};
            
            std::vector<double> genotype_log_posteriors(genotypes.size());
            
            std::transform(std::cbegin(genotypes), std::cend(genotypes), std::begin(genotype_log_posteriors),
                           [&] (const auto& genotype) {
                               std::vector<Haplotype> haps {std::cbegin(genotype), std::cend(genotype)};
                               return std::log(haplotype_model.evaluate(haps))
                                            + likelihood_model.log_likelihood(sample, genotype);
                           });
            
            const auto norm = Maths::log_sum_exp(genotype_log_posteriors);
            
            for (auto& p : genotype_log_posteriors) p = std::exp(p -= norm);
            
            auto haplotype_posteriors = calculate_haplotype_posteriors(haplotypes, genotypes,
                                                                       genotype_log_posteriors);
            
            ProbabilityMatrix<Genotype<Haplotype>> gps {
                std::make_move_iterator(std::begin(genotypes)),
                std::make_move_iterator(std::end(genotypes))
            };
            
            insert_sample(sample, genotype_log_posteriors, gps);
            
            return Latents {std::move(gps), std::move(haplotype_posteriors)};
        }
    } // namesapce GenotypeModel
} // namespace Octopus
