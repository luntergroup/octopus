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
#include "equitable.h"

//using Genotype = std::vector<Haplotype>;

class Genotype : public Equitable<Genotype>
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
    std::vector<Haplotype> get_unique_haplotypes() const;
    //const Haplotype& get_first_alternate_haplotype(const Haplotype& a_haplotype) const;
    
    friend bool operator==(const Genotype& lhs, const Genotype& rhs);
private:
    std::vector<Haplotype> the_haplotypes_;
};

unsigned num_genotypes(unsigned num_haplotypes, unsigned ploidy);

// Assumes the input haplotypes are all unique
std::vector<Genotype> get_all_genotypes(const std::vector<Haplotype>& haplotypes, unsigned ploidy);

inline bool operator==(const Genotype& lhs, const Genotype& rhs)
{
    if (lhs.ploidy() != rhs.ploidy()) return false;
    
    for (const auto& lhs_haplotype : lhs.the_haplotypes_) {
        if (lhs.num_occurences(lhs_haplotype) != rhs.num_occurences(lhs_haplotype)) return false;
    }
    
    return true;
}

namespace std {
    template <> struct hash<Genotype>
    {
        size_t operator()(const Genotype& g) const
        {
            return hash<Haplotype>()(g.at(0)); //TODO: improve this
        }
    };
}

#endif /* defined(__Octopus__genotype__) */
