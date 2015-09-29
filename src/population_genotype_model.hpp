//
//  population_genotype_model.hpp
//  Octopus
//
//  Created by Daniel Cooke on 26/08/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__population_genotype_model__
#define __Octopus__population_genotype_model__

#include "genotype_model.hpp"

#include <unordered_map>

namespace Octopus
{
    class PopulationGenotypeModel : public GenotypeModel
    {
    public:
        PopulationGenotypeModel(unsigned num_samples, unsigned sample_ploidy,
                                unsigned max_em_iterations = 100, double em_epsilon = 0.001);
        
    private:
        using HaplotypeFrequencies = std::unordered_map<Haplotype, double>;
        
        const unsigned max_em_iterations_;
        const double em_epsilon_;
        
        const unsigned num_samples_;
        const unsigned sample_ploidy_;
        
        GenotypeProbabilities do_evaluate(const std::vector<Haplotype>& haplotypes, const ReadMap& reads) override;
    };
} // namespace Octopus

#endif /* defined(__Octopus__population_genotype_model__) */
