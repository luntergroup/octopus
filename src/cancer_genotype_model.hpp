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

#include "genotype_model.hpp"
#include "mappable_map.hpp"
#include "cancer_genotype.hpp"
#include "reference_genome.hpp"
#include "haplotype_prior_model.hpp"
#include "single_read_model.hpp"

class AlignedRead;
class Haplotype;

namespace Octopus
{
    namespace GenotypeModel
    {
    class Cancer
    {
    public:
        using GenotypeProbabilities          = std::unordered_map<CancerGenotype<Haplotype>, double>;
        using GenotypeWeightResponsibilities = std::unordered_map<SampleIdType, std::vector<std::array<double, 3>>>;
        using SampleGenotypeWeightsPriors    = std::array<double, 3>;
        using SampleGenotypeWeights          = std::array<double, 3>;
        using GenotypeWeightsPriors          = std::unordered_map<SampleIdType, SampleGenotypeWeightsPriors>;
        using GenotypeWeights                = std::unordered_map<SampleIdType, SampleGenotypeWeights>;
        
        struct Latents
        {
            Latents() = default;
            template <typename G, typename W>
            Latents(G&& genotype_posteriors, W&& genotype_weights)
            :
            genotype_posteriors {std::forward<G>(genotype_posteriors)},
            genotype_weights {std::forward<W>(genotype_weights)}
            {}
            
            GenotypeProbabilities genotype_posteriors;
            GenotypeWeights genotype_weights;
        };
        
        Cancer(SampleIdType normal_sample, unsigned max_em_iterations = 100, double em_epsilon = 0.001);
        
        Latents evaluate(const std::vector<Haplotype>& haplotypes, const ReadMap& reads, ReferenceGenome& reference);
        
    private:
        HaplotypePriorModel haplotype_prior_model_;
        SingleReadModel read_model_;
        
        unsigned max_em_iterations_;
        double em_epsilon_;
        
        SampleIdType normal_sample_;
    };
    
    } // namespace GenotypeModel
} // namespace Octopus

#endif /* defined(__Octopus__cancer_genotype_model__) */
