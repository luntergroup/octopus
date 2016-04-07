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

#include "common.hpp"
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
        using HaplotypePriorMap  = std::unordered_map<HaplotypeReference, double>;
        
        struct Latents
        {
            using HaplotypeFrequencyMap   = std::unordered_map<HaplotypeReference, double>;
            using GenotypeProbabilityMap  = ProbabilityMatrix<Genotype<Haplotype>>;
            using HaplotypeProbabilityMap = std::unordered_map<HaplotypeReference, double>;
            
            Latents() = default;
            template <typename F, typename G, typename P>
            Latents(F&& haplotype_frequencies, G&& genotype_posteriors, P&& haplotype_posteriors);
            ~Latents() = default;
            
            HaplotypeFrequencyMap haplotype_frequencies;
            GenotypeProbabilityMap genotype_posteriors;
            HaplotypeProbabilityMap haplotype_posteriors;
        };
        
        explicit Population(unsigned ploidy, unsigned max_em_iterations = 100, double em_epsilon = 0.001);
        
        Latents infer_latents(const std::vector<SampleIdType>& samples,
                              const std::vector<Haplotype>& haplotypes,
                              const HaplotypePriorMap& haplotype_priors,
                              const HaplotypeLikelihoodCache& haplotype_likelihoods) const;
        
    private:
        const unsigned ploidy_;
        const unsigned max_em_iterations_;
        const double em_epsilon_;
    };
    
    template <typename F, typename G, typename P>
    Population::Latents::Latents(F&& haplotype_frequencies, G&& genotype_posteriors,
                                 P&& haplotype_posteriors)
    :
    haplotype_frequencies {std::forward<F>(haplotype_frequencies)},
    genotype_posteriors {std::forward<G>(genotype_posteriors)},
    haplotype_posteriors {std::forward<P>(haplotype_posteriors)}
    {}
        
    } // namespace GenotypeModel
} // namespace Octopus

#endif /* defined(__Octopus__population_genotype_model__) */
