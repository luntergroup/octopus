//
//  genotype.h
//  Octopus
//
//  Created by Daniel Cooke on 24/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__genotype__
#define __Octopus__genotype__

#include <vector>

#include "haplotype.h"

using Genotype = std::vector<Haplotype>;

bool is_homozygous(const Genotype& a_genotype);

bool is_heterozygous(const Genotype& a_genotype);

unsigned num_occurences(const Haplotype& a_haplotype, const Genotype& a_genotype);

unsigned ploidy(const Genotype& a_genotype);

const Haplotype& get_first_alternate_haplotype(const Haplotype& a_haplotype, const Genotype& a_genotype);

// Assumes the input haplotypes are all unique
std::vector<Genotype> get_all_genotypes(unsigned ploidy, const std::vector<Haplotype>& haplotypes);

#endif /* defined(__Octopus__genotype__) */
