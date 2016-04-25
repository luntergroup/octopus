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
                using GenotypeProbabilityVector = std::vector<double>;
                
                Latents() = default;
                Latents(GenotypeProbabilityVector&& genotype_probabilities);
                ~Latents() = default;
                
                GenotypeProbabilityVector genotype_probabilities;
            };
            
            struct InferredLatents
            {
                InferredLatents() = default;
                InferredLatents(Latents&& posteriors, double log_evidence);
                
                Latents posteriors;
                double log_evidence;
            };
            
            Individual() = delete;
            
            explicit Individual(unsigned ploidy, const CoalescentModel& genotype_prior_model);
            
            ~Individual() = default;
            
            Individual(const Individual&)            = delete;
            Individual& operator=(const Individual&) = delete;
            Individual(Individual&&)                 = delete;
            Individual& operator=(Individual&&)      = delete;
            
            InferredLatents infer_latents(const SampleIdType& sample,
                                          const std::vector<Genotype<Haplotype>>& genotypes,
                                          const HaplotypeLikelihoodCache& haplotype_likelihoods) const;
            
        private:
            unsigned ploidy_;
            std::reference_wrapper<const CoalescentModel> genotype_prior_model_;
        };
    } // namesapce GenotypeModel
} // namespace Octopus

#endif /* individual_genotype_model_hpp */
