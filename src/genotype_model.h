//
//  genotype_model.h
//  Octopus
//
//  Created by Daniel Cooke on 22/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__genotype_model__
#define __Octopus__genotype_model__

#include <vector>
#include <cstddef>

#include "haplotype.h"

using std::size_t;

class AlignedRead;

class GenotypeModel
{
public:
    using SampleProbabilities = std::vector<double>;
    struct HaplotypePosteriors
    {
        double population_probability;
        SampleProbabilities sample_probabilities;
    };
    using Haplotypes  = std::vector<Haplotype>;
    using SampleReads = std::vector<std::vector<AlignedRead>>;
    
    GenotypeModel() = delete;
    explicit GenotypeModel(unsigned ploidy, size_t max_num_haplotypes);
    ~GenotypeModel() = default;
    
    GenotypeModel(const GenotypeModel&)            = default;
    GenotypeModel& operator=(const GenotypeModel&) = default;
    GenotypeModel(GenotypeModel&&)                 = default;
    GenotypeModel& operator=(GenotypeModel&&)      = default;
    
    std::vector<HaplotypePosteriors> get_haplotype_probabilities(const Haplotypes& the_haplotypes,
                                                                 const SampleReads& the_reads);
private:
    using Genotype               = std::vector<Haplotype>;
    using HaplotypeProbabilities = std::vector<double>;
    
    unsigned ploidy_;
    size_t max_num_haplotypes_;
    
    double genotype_probability(const Genotype& the_genotype,
                                const HaplotypeProbabilities& the_haplotype_probabilities) const;
    double read_probability(const AlignedRead& a_read, const Genotype& the_genotype) const;
    double reads_liklihood(const std::vector<AlignedRead>& individual_reads,
                           const Genotype& the_genotype) const;
    double genotype_likilihood(const std::vector<AlignedRead>& individual_reads,
                               const Genotype& the_genotype,
                               const HaplotypeProbabilities& the_haplotype_probabilities) const;
};

#endif /* defined(__Octopus__genotype_model__) */
