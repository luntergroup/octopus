//
//  individual_genotype_model.hpp
//  Octopus
//
//  Created by Daniel Cooke on 01/04/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#ifndef individual_genotype_model_hpp
#define individual_genotype_model_hpp

#include <vector>
#include <unordered_map>
#include <functional>

#include "common.hpp"
#include "haplotype.hpp"
#include "coalescent_model.hpp"
#include "haplotype_likelihood_cache.hpp"
#include "genotype.hpp"
#include "probability_matrix.hpp"

namespace Octopus
{
    namespace GenotypeModel
    {
        class Individual
        {
        public:
            struct Latents
            {
                using GenotypeProbabilityMap = ProbabilityMatrix<Genotype<Haplotype>>;
                
                Latents() = default;
                Latents(GenotypeProbabilityMap&& genotype_posteriors);
                ~Latents() = default;
                
                GenotypeProbabilityMap genotype_posteriors;
            };
            
            Individual()  = delete;
            explicit Individual(unsigned ploidy);
            ~Individual() = default;
            
            Individual(const Individual&)            = delete;
            Individual& operator=(const Individual&) = delete;
            Individual(Individual&&)                 = delete;
            Individual& operator=(Individual&&)      = delete;
            
            Latents
            infer_latents(const SampleIdType& sample,
                          std::vector<Genotype<Haplotype>> candidate_genotypes,
                          const CoalescentModel& genotype_prior_model,
                          const HaplotypeLikelihoodCache& haplotype_likelihoods) const;
            
            Latents
            infer_latents(const SampleIdType& sample,
                          const std::vector<Haplotype>& candidate_haplotypes,
                          const CoalescentModel& genotype_prior_model,
                          const HaplotypeLikelihoodCache& haplotype_likelihoods) const;
            
            double log_evidence(const SampleIdType& sample,
                                const HaplotypeLikelihoodCache& haplotype_likelihoods,
                                const Latents& latents) const;
            
        private:
            unsigned ploidy_;
        };
    } // namesapce GenotypeModel
} // namespace Octopus

#endif /* individual_genotype_model_hpp */
