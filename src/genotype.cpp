//
//  genotype.cpp
//  Octopus
//
//  Created by Daniel Cooke on 24/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "genotype.h"

//#include <algorithm> // std::inplace_merge, std::all_of, std::binary_search, std::equal_range, std::unique_copy
#include <boost/math/special_functions/binomial.hpp>

const Haplotype& Genotype::at(unsigned n) const
{
    return the_haplotypes_.at(n);
}

void Genotype::emplace(const Haplotype& haplotype)
{
    the_haplotypes_.emplace_back(haplotype);
    std::inplace_merge(std::begin(the_haplotypes_), std::prev(std::end(the_haplotypes_)), std::end(the_haplotypes_));
}

void Genotype::emplace(Haplotype&& haplotype)
{
    the_haplotypes_.emplace_back(std::move(haplotype));
    std::inplace_merge(std::begin(the_haplotypes_), std::prev(std::end(the_haplotypes_)), std::end(the_haplotypes_));
}

bool Genotype::is_homozygous() const
{
    const auto& first_haplotype = the_haplotypes_.front();
    return std::all_of(std::next(std::cbegin(the_haplotypes_)), std::cend(the_haplotypes_),
                       [&first_haplotype] (const auto& a_haplotype) {
                           return first_haplotype == a_haplotype;
                       });
}

bool Genotype::contains(const Haplotype& a_haplotype) const
{
    return std::binary_search(std::cbegin(the_haplotypes_), std::cend(the_haplotypes_), a_haplotype);
}

unsigned Genotype::num_occurences(const Haplotype& a_haplotype) const
{
    auto equal_range = std::equal_range(std::cbegin(the_haplotypes_), std::cend(the_haplotypes_), a_haplotype);
    return static_cast<unsigned>(std::distance(equal_range.first, equal_range.second));
}

unsigned Genotype::ploidy() const noexcept
{
    return static_cast<unsigned>(the_haplotypes_.size());
}

std::vector<Haplotype> Genotype::get_unique_haplotypes() const
{
    std::vector<Haplotype> result {};
    result.reserve(ploidy());
    
    std::unique_copy(std::cbegin(the_haplotypes_), std::cend(the_haplotypes_), std::back_inserter(result));
    
    return result;
}

unsigned num_genotypes(unsigned num_haplotypes, unsigned ploidy)
{
    return static_cast<unsigned>(boost::math::binomial_coefficient<double>(num_haplotypes + ploidy - 1,
                                                                           num_haplotypes - 1));
}

Genotype get_genotype_from_haplotype_indicies(const std::vector<Haplotype>& haplotypes,
                                              const std::vector<unsigned>& haplotype_indicies)
{
    Genotype result {};
    
    for (auto i : haplotype_indicies) {
        result.emplace(haplotypes.at(i));
    }
    
    return result;
}

// This is a surprisingly complex algorithm. Essentially we are enumerating the terms in a
// polynomial of the form (x_1 + ... + x_{num_haplotypes})^ploidy. It's easier if we first map
// haplotypes to integers such that h_1 -> 0, ..., h_num_haplotypes -> num_haplotypes - 1.
// The order of the haplotypes is not important, only that each haplotype maps to a unique integer,
// and is within the set {0, num_haplotypes - 1}.
// We then start with the homozygous genotype {0, ..., 0} and increment this each time
// we exit the loop. We are then essentially enumerating integers in base_{num_haplotypes}.
// However there is a caveat, usually when we increment an integer and need to carry, the bits to the
// right of the carried 1 are all set to zero (e.g. in base 3, 0122 + 1 = 0200). In our case this
// doesn't work as we will already have the resulting genotype (0002 in the example). The solution
// is to 'miss out' all the genotypes we already have by setting the bits to the right of the carry
// to the same value as the left most significant bit (e.g. 0122 + 1 = 0222). This is guaranteed
// to be a novel genotype. We terminate when the leftmost bit is equal to the num_haplotypes - 1
// which is necessarily the genotype {num_haplotypes - 1, ..., num_haplotypes - 1}.
std::vector<Genotype> get_all_genotypes(const std::vector<Haplotype>& haplotypes, unsigned ploidy)
{
    std::vector<Genotype> result {};
    
    unsigned num_haplotypes {static_cast<unsigned>(haplotypes.size())};
    
    std::vector<unsigned> haplotype_indicies(ploidy, 0);
    
    unsigned i {};
    while (true) {
        if (haplotype_indicies[i] == num_haplotypes) {
            do {
                ++i;
            } while (haplotype_indicies[i] == num_haplotypes - 1);
            
            if (i == ploidy) return result;
            
            ++haplotype_indicies[i];
            
            for (unsigned j {0}; j <= i; ++j) {
                haplotype_indicies[j] = haplotype_indicies[i];
            }
            
            i = 0;
        }
        
        result.push_back(get_genotype_from_haplotype_indicies(haplotypes, haplotype_indicies));
        ++haplotype_indicies[i];
    }
}

std::unordered_map<Haplotype, unsigned> get_haplotype_occurence_map(const Genotype& a_genotype)
{
    std::unordered_map<Haplotype, unsigned> result {};
    
    for (unsigned i {}; i < a_genotype.ploidy(); ++i) {
        ++result[a_genotype.at(i)];
    }
    
    return result;
}
