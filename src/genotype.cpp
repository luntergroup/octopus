//
//  genotype.cpp
//  Octopus
//
//  Created by Daniel Cooke on 24/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "genotype.h"

#include <algorithm> // std::all_of, std::count, std::find
#include <iterator>  // std::cbegin etc
#include <boost/math/special_functions/binomial.hpp>

#include "utils.h"

const Haplotype& Genotype::at(unsigned n) const
{
    return the_haplotypes_.at(n);
}

void Genotype::emplace(const Haplotype& haplotype)
{
    the_haplotypes_.emplace_back(haplotype);
}

void Genotype::emplace(Haplotype&& haplotype)
{
    the_haplotypes_.emplace_back(std::move(haplotype));
}

bool Genotype::is_homozygous() const
{
    const auto& first_haplotype = the_haplotypes_.front();
    return std::all_of(std::next(std::cbegin(the_haplotypes_)), std::cend(the_haplotypes_),
                       [&first_haplotype] (const auto& a_haplotype) {
                           return first_haplotype == a_haplotype;
                       });
}

unsigned Genotype::num_occurences(const Haplotype& a_haplotype) const
{
    return static_cast<unsigned>(std::count(std::cbegin(the_haplotypes_), std::cend(the_haplotypes_),
                                            a_haplotype));
}

unsigned Genotype::ploidy() const noexcept
{
    return static_cast<unsigned>(the_haplotypes_.size());
}

std::vector<Haplotype> Genotype::get_unique_haplotypes() const
{
    std::vector<Haplotype> result {};
    
    for (const auto& haplotype : the_haplotypes_) {
        if (std::find(result.cbegin(), result.cend(), haplotype) == result.cend()) {
            result.push_back(haplotype);
        }
    }
    
    return result;
}

//const Haplotype& Genotype::get_first_alternate_haplotype(const Haplotype& a_haplotype) const
//{
//    return *std::find_if_not(std::cbegin(the_haplotypes_), std::cend(the_haplotypes_), a_haplotype);
//}

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

std::vector<Genotype> get_all_genotypes(const std::vector<Haplotype>& haplotypes, unsigned ploidy)
{
    std::vector<Genotype> result {};
    
    unsigned num_haplotypes {static_cast<unsigned>(haplotypes.size())};
    
    std::vector<unsigned> haplotype_indicies(ploidy, 0);
    
    unsigned i {0};
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
