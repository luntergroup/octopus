//
//  population_genotype_model.h
//  Octopus
//
//  Created by Daniel Cooke on 26/08/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__population_genotype_model__
#define __Octopus__population_genotype_model__

#include "genotype_model.h"

#include <unordered_map>

namespace Octopus
{
    class PopulationGenotypeModel : public GenotypeModel
    {
    public:
        PopulationGenotypeModel(unsigned num_samples, unsigned sample_ploidy);
        
    private:
        using HaplotypeFrequencies = std::unordered_map<Haplotype, double>;
        
        unsigned num_samples_;
        unsigned sample_ploidy_;
        
        GenotypeProbabilities do_evaluate(const std::vector<Haplotype>& haplotypes, const ReadMap& reads) override;
    };
} // end namespace Octopus

#endif /* defined(__Octopus__population_genotype_model__) */
