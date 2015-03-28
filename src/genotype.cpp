//
//  genotype.cpp
//  Octopus
//
//  Created by Daniel Cooke on 24/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "genotype.h"

#include <algorithm> // std::all_of, std::count
#include <iterator>  // std::cbegin etc

#include "utils.h"

bool is_homozygous(const Genotype& a_genotype)
{
    const auto& first_haplotype = a_genotype.front();
    return std::all_of(std::next(std::cbegin(a_genotype)), std::cend(a_genotype),
                       [first_haplotype] (const auto& a_haplotype) {
                           first_haplotype == a_haplotype;
                       });
}

bool is_heterozygous(const Genotype& a_genotype)
{
    return !is_homozygous(a_genotype);
}

unsigned num_occurences(const Haplotype& a_haplotype, const Genotype& a_genotype)
{
    return static_cast<unsigned>(std::count(std::cbegin(a_genotype), std::cend(a_genotype), a_haplotype));
}

unsigned ploidy(const Genotype& a_genotype)
{
    return static_cast<unsigned>(a_genotype.size());
}

const Haplotype& get_first_alternate_haplotype(const Haplotype& a_haplotype, const Genotype& a_genotype)
{
    return *std::find_if_not(std::cbegin(a_genotype), std::cend(a_genotype), a_haplotype);
}

std::vector<Genotype> get_all_genotypes(unsigned ploidy, const std::vector<Haplotype>& haplotypes)
{
    
    
    return std::vector<Genotype> {};
}
