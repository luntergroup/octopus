//
//  population_genotype_model.hpp
//  Octopus
//
//  Created by Daniel Cooke on 26/08/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__population_genotype_model__
#define __Octopus__population_genotype_model__

#include <string>
#include <vector>
#include <unordered_map>
#include <functional>
#include <utility>

#include "genotype.hpp"
#include "reference_genome.hpp"
#include "haplotype_likelihood_cache.hpp"
#include "probability_matrix.hpp"

class AlignedRead;
class Haplotype;

namespace Octopus
{
    namespace GenotypeModel
    {
    class Population
    {
    public:
        using HaplotypeReference = std::reference_wrapper<const Haplotype>;
        using HaplotypePrioMap   = std::unordered_map<HaplotypeReference, double>;
        
        struct Latents
        {
            using HaplotypeFrequencyMap  = std::unordered_map<HaplotypeReference, double>;
            using GenotypeProbabilityMap = ProbabilityMatrix<Genotype<Haplotype>>;
            
            Latents() = default;
            template <typename G, typename F>
            Latents(G&& genotype_posteriors, F&& haplotype_frequencies)
            :
            genotype_posteriors {std::forward<G>(genotype_posteriors)},
            haplotype_frequencies {std::forward<F>(haplotype_frequencies)}
            {}
            
            GenotypeProbabilityMap genotype_posteriors;
            HaplotypeFrequencyMap haplotype_frequencies;
        };
        
        explicit Population(unsigned ploidy, unsigned max_em_iterations = 100, double em_epsilon = 0.001);
        
        Latents infer_latents(const std::vector<Haplotype>& haplotypes,
                              const HaplotypePrioMap& haplotype_priors,
                              HaplotypeLikelihoodCache& haplotype_likelihoods,
                              const ReadMap& reads);
        
    private:
        const unsigned ploidy_;
        const unsigned max_em_iterations_;
        const double em_epsilon_;
    };
    
    } // namespace GenotypeModel
} // namespace Octopus

#endif /* defined(__Octopus__population_genotype_model__) */
