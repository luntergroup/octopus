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
#include <unordered_map>
#include <ostream>
#include <iterator>
#include <boost/functional/hash.hpp> // boost::hash_combine

#include "haplotype.h"
#include "equitable.h"

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
    bool contains(const Haplotype& a_haplotype) const;
    unsigned num_occurences(const Haplotype& a_haplotype) const;
    bool is_homozygous() const;
    std::vector<Haplotype> get_unique_haplotypes() const;
    
    friend bool operator==(const Genotype& lhs, const Genotype& rhs);
private:
    std::vector<Haplotype> the_haplotypes_;
};

unsigned num_genotypes(unsigned num_haplotypes, unsigned ploidy);

// Assumes the input haplotypes are all unique
std::vector<Genotype> get_all_genotypes(const std::vector<Haplotype>& haplotypes, unsigned ploidy);

std::unordered_map<Haplotype, unsigned> get_haplotype_occurence_map(const Genotype& a_genotype);

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
            size_t seed {};
            
            for (unsigned i {}; i < g.ploidy(); ++i) {
                boost::hash_combine(seed, hash<Haplotype>()(g.at(i)));
            }
            
            return seed;
        }
    };
}

inline std::ostream& operator<<(std::ostream& os, const Genotype& a_genotype)
{
    auto haplotype_occurences = get_haplotype_occurence_map(a_genotype);
    std::vector<std::pair<Haplotype, unsigned>> p {haplotype_occurences.begin(), haplotype_occurences.end()};
    for (unsigned i {}; i < p.size() - 1; ++i) {
        os << p[i].first << "(" << p[i].second << "),";
    }
    os << p.back().first << "(" << p.back().second << ")";
    return os;
}

#endif /* defined(__Octopus__genotype__) */
