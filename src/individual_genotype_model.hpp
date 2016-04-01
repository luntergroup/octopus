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
                using HaplotypeReference      = std::reference_wrapper<const Haplotype>;
                using GenotypeProbabilityMap  = ProbabilityMatrix<Genotype<Haplotype>>;
                using HaplotypeProbabilityMap = std::unordered_map<HaplotypeReference, double>;
                
                Latents() = default;
                
                template <typename G, typename P>
                Latents(G&& genotype_posteriors, P&& haplotype_posteriors);
                
                ~Latents() = default;
                
                GenotypeProbabilityMap genotype_posteriors;
                HaplotypeProbabilityMap haplotype_posteriors;
            };
            
            Individual()  = delete;
            explicit Individual(unsigned ploidy);
            ~Individual() = default;
            
            Individual(const Individual&)            = delete;
            Individual& operator=(const Individual&) = delete;
            Individual(Individual&&)                 = delete;
            Individual& operator=(Individual&&)      = delete;
            
            Latents infer_latents(const SampleIdType& sample,
                                  const std::vector<Haplotype>& candidate_haplotypes,
                                  const CoalescentModel& haplotype_model,
                                  const HaplotypeLikelihoodCache& haplotype_likelihoods) const;
        private:
            unsigned ploidy_;
        };
        
        template <typename G, typename P>
        Individual::Latents::Latents(G&& genotype_posteriors, P&& haplotype_posteriors)
        :
        genotype_posteriors {std::forward<G>(genotype_posteriors)},
        haplotype_posteriors {std::forward<P>(haplotype_posteriors)}
        {}
    } // namesapce GenotypeModel
} // namespace Octopus

#endif /* individual_genotype_model_hpp */
