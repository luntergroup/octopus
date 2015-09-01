//
//  population_genotype_model.cpp
//  Octopus
//
//  Created by Daniel Cooke on 26/08/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "population_genotype_model.h"

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
    
    double log_hardy_weinberg(const Genotype<Haplotype>& genotype,
                              const std::unordered_map<Haplotype, double>& haplotype_frequencies)
    {
        auto unique_haplotypes = genotype.get_unique();
        
        std::vector<unsigned> occurences {};
        occurences.reserve(unique_haplotypes.size());
        
        double r {};
        
        for (const auto& haplotype : unique_haplotypes) {
            auto num_occurences = genotype.num_occurences(haplotype);
            occurences.push_back(num_occurences);
            r += num_occurences * haplotype_frequencies.at(haplotype);
        }
        
        auto c = log_multinomial_coefficient<double>(occurences.cbegin(), occurences.cend());
        
        return c * r;
    }
    
} // end namespace Octopus

