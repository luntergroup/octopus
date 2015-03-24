//
//  genotype_model.cpp
//  Octopus
//
//  Created by Daniel Cooke on 22/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "genotype_model.h"

#include "maths.h"

GenotypeModel::GenotypeModel(unsigned ploidy, size_t max_num_haplotypes)
:
ploidy_ {ploidy},
max_num_haplotypes_ {max_num_haplotypes}
{}

std::vector<GenotypeModel::HaplotypePosteriors>
GenotypeModel::get_haplotype_probabilities(const Haplotypes& the_haplotypes,
                                           const SampleReads& the_reads)
{
    return std::vector<HaplotypePosteriors> {};
}

double GenotypeModel::genotype_probability(const Genotype& the_genotype,
                                           const HaplotypeProbabilities& the_haplotype_probabilities) const
{
    return 0;
}

double GenotypeModel::read_probability(const AlignedRead& a_read, const Genotype& the_genotype) const
{
    return 0;
}

double GenotypeModel::reads_liklihood(const std::vector<AlignedRead>& individual_reads,
                                      const Genotype& the_genotype) const
{
    double result {1};
    
    for (const auto& read : individual_reads) {
        result *= read_probability(read, the_genotype);
    }
    
    return result;
}

double GenotypeModel::genotype_likilihood(const std::vector<AlignedRead>& individual_reads,
                                          const Genotype& the_genotype,
                                          const HaplotypeProbabilities& the_haplotype_probabilities) const
{
    return genotype_probability(the_genotype, the_haplotype_probabilities) *
                reads_liklihood(individual_reads, the_genotype);
}
