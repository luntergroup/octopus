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

#include "genotype_model.hpp"
#include "mappable_map.hpp"
#include "genotype.hpp"

class AlignedRead;
class Haplotype;

namespace Octopus
{
    namespace GenotypeModel
    {
    class Population
    {
    public:
        using SampleGenotypeProbabilities = std::unordered_map<Genotype<Haplotype>, double>;
        using GenotypeProbabilities       = std::unordered_map<SampleIdType, SampleGenotypeProbabilities>;
        
        struct Latents
        {
            GenotypeProbabilities genotype_posteriors;
            HaplotypeFrequencies haplotype_frequencies;
        };
        
        Population(unsigned ploidy, unsigned max_em_iterations = 100, double em_epsilon = 0.001);
        
        Latents evaluate(const std::vector<Haplotype>& haplotypes, const ReadMap& reads);
        
    private:
        const unsigned ploidy_;
        const unsigned max_em_iterations_;
        const double em_epsilon_;
    };
    
    } // namespace GenotypeModel
} // namespace Octopus

#endif /* defined(__Octopus__population_genotype_model__) */
