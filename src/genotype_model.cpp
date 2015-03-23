//
//  genotype_model.cpp
//  Octopus
//
//  Created by Daniel Cooke on 22/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "genotype_model.h"

GenotypeModel::GenotypeModel(unsigned ploidy, size_t max_num_haplotypes)
:
ploidy_ {ploidy},
max_num_haplotypes_ {max_num_haplotypes}
{}

std::vector<GenotypeModel::HaplotypeProbabilities>
GenotypeModel::get_haplotype_probabilities(const Haplotypes& the_haplotypes,
                                           const SampleReads& the_reads)
{
    return std::vector<HaplotypeProbabilities> {};
}
