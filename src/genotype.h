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
#include <iterator> // std::cbegin etc
#include <initializer_list>
#include <boost/functional/hash.hpp> // boost::hash_range

#include "haplotype.h"
#include "equitable.h"

class Genotype : public Equitable<Genotype>
{
public:
    using HaplotypeIterator = std::vector<Haplotype>::const_iterator;
    
    Genotype() = default;
    Genotype(unsigned ploidy);
    Genotype(std::initializer_list<Haplotype> haplotypes);
    ~Genotype() = default;
    
    Genotype(const Genotype&)            = default;
    Genotype& operator=(const Genotype&) = default;
    Genotype(Genotype&&)                 = default;
    Genotype& operator=(Genotype&&)      = default;
    
    const Haplotype& at(unsigned n) const;
    void emplace(const Haplotype& haplotype);
    void emplace(Haplotype&& haplotype);
    
    HaplotypeIterator begin() const;
    HaplotypeIterator end() const;
    HaplotypeIterator cbegin() const;
    HaplotypeIterator cend() const;
    
    unsigned ploidy() const noexcept;
    bool contains(const Haplotype& a_haplotype) const;
    unsigned num_occurences(const Haplotype& a_haplotype) const;
    bool is_homozygous() const;
    std::vector<Haplotype> get_unique_haplotypes() const;
    
    friend bool operator==(const Genotype& lhs, const Genotype& rhs);
    
private:
    std::vector<Haplotype> the_haplotypes_;
};

bool operator==(const Genotype& lhs, const Genotype& rhs);

namespace std {
    template <> struct hash<Genotype>
    {
        size_t operator()(const Genotype& genotype) const
        {
            return boost::hash_range(std::cbegin(genotype), std::cend(genotype));
        }
    };
}

unsigned num_genotypes(unsigned num_haplotypes, unsigned ploidy);

// Assumes the input haplotypes are all unique
std::vector<Genotype> get_all_genotypes(const std::vector<Haplotype>& haplotypes, unsigned ploidy);

std::unordered_map<Haplotype, unsigned> get_haplotype_occurence_map(const Genotype& a_genotype);

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
