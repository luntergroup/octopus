//
//  population_genotype_model.cpp
//  Octopus
//
//  Created by Daniel Cooke on 26/08/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "population_genotype_model.h"

#include <unordered_map>

#include "read_model.h"
#include "maths.h"

namespace Octopus
{
    // public methods
    
    PopulationGenotypeModel::PopulationGenotypeModel(unsigned num_samples, unsigned sample_ploidy)
    :
    num_samples_ {num_samples},
    sample_ploidy_ {sample_ploidy}
    {}
    
    // private methods
    
    PopulationGenotypeModel::GenotypeProbabilities
    PopulationGenotypeModel::do_evaluate(const std::vector<Haplotype>& haplotypes, const ReadMap& reads)
    {
        GenotypeProbabilities result {};
        
        return result;
    }
    
    // non member methods
    
    double hardy_weinberg_log_probability(const Genotype<Haplotype>& genotype,
                                          const std::unordered_map<Haplotype, double>& haplotype_frequencies)
    {
        return 0;
    }
    
} // end namespace Octopus

