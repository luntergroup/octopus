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

//using Genotype = std::vector<Haplotype>;

class Genotype
{
public:
    Genotype()  = default;
    ~Genotype() = default;
    
    Genotype(const Genotype&)            = default;
    Genotype& operator=(const Genotype&) = default;
    Genotype(Genotype&&)                 = default;
    Genotype& operator=(Genotype&&)      = default;
    
    const Haplotype& at(unsigned n) const;
    void emplace(const Haplotype& haplotype);
    void emplace(Haplotype&& haplotype);
    
    unsigned ploidy() const noexcept;
    unsigned num_occurences(const Haplotype& a_haplotype) const;
    bool is_homozygous() const;
    //const Haplotype& get_first_alternate_haplotype(const Haplotype& a_haplotype) const;

private:
    std::vector<Haplotype> the_haplotypes_;
};

unsigned num_genotypes(unsigned num_haplotypes, unsigned ploidy);

// Assumes the input haplotypes are all unique
std::vector<Genotype> get_all_genotypes(const std::vector<Haplotype>& haplotypes, unsigned ploidy);

#endif /* defined(__Octopus__genotype__) */
