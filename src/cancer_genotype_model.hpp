//
//  cancer_genotype_model.hpp
//  Octopus
//
//  Created by Daniel Cooke on 26/08/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__cancer_genotype_model__
#define __Octopus__cancer_genotype_model__

#include <string>
#include <vector>
#include <unordered_map>

#include "common.hpp"
#include "cancer_genotype.hpp"
#include "reference_genome.hpp"
#include "map_utils.hpp"
#include "haplotype_likelihood_cache.hpp"

class AlignedRead;
class Haplotype;

namespace Octopus
{
    namespace GenotypeModel
    {
    class Cancer
    {
    public:
        using HaplotypeReference = std::reference_wrapper<const Haplotype>;
        using HaplotypePrioMap   = std::unordered_map<HaplotypeReference, double>;
        
        using GenotypeProbabilityMap       = std::unordered_map<CancerGenotype<Haplotype>, double>;
        using SampleGenotypeMixturesPriors = std::array<double, 3>;
        using SampleGenotypeMixtures       = std::array<double, 3>;
        using GenotypeMixturesPriors       = std::unordered_map<SampleIdType, SampleGenotypeMixturesPriors>;
        using GenotypeMixtures             = std::unordered_map<SampleIdType, SampleGenotypeMixtures>;
        
        struct Latents
        {
            Latents() = default;
            template <typename G, typename W>
            Latents(G&& genotype_posteriors, W&& genotype_mixtures)
            :
            genotype_posteriors {std::forward<G>(genotype_posteriors)},
            genotype_mixtures {std::forward<W>(genotype_mixtures)}
            {}
            
            GenotypeProbabilityMap genotype_posteriors;
            GenotypeMixtures genotype_mixtures;
        };
        
        explicit Cancer(SampleIdType normal_sample, unsigned max_em_iterations = 100,
                        double em_epsilon = 0.001);
        
        Latents infer_latents(const std::vector<Haplotype>& haplotypes,
                              const HaplotypePrioMap& haplotype_priors,
                              HaplotypeLikelihoodCache& haplotype_likelihoods,
                              const ReadMap& reads);
        
    private:
        unsigned max_em_iterations_;
        double em_epsilon_;
        
        SampleIdType normal_sample_;
    };
    
    } // namespace GenotypeModel
} // namespace Octopus

#endif /* defined(__Octopus__cancer_genotype_model__) */
