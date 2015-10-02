//
//  genotype.cpp
//  Octopus
//
//  Created by Daniel Cooke on 02/10/2015.
//  Copyright © 2015 Oxford University. All rights reserved.
//

#include "genotype.hpp"

#include <functional>  // std::not_equal_to
#include <boost/math/special_functions/binomial.hpp>

Genotype<Allele>::Genotype(unsigned ploidy)
:
alleles_ {}
{
    alleles_.reserve(ploidy);
}

Genotype<Allele>::Genotype(unsigned ploidy, const Allele& init)
:
alleles_ {ploidy, init}
{}

Genotype<Allele>::Genotype(std::initializer_list<Allele> alleles)
:
alleles_ {alleles}
{}

const Allele& Genotype<Allele>::at(unsigned n) const
{
    return alleles_.at(n);
}

const Allele& Genotype<Allele>::operator[](unsigned n) const
{
    return alleles_[n];
}

typename Genotype<Allele>::Iterator Genotype<Allele>::begin() const noexcept
{
    return alleles_.begin();
}

typename Genotype<Allele>::Iterator Genotype<Allele>::end() const noexcept
{
    return alleles_.end();
}

typename Genotype<Allele>::Iterator Genotype<Allele>::cbegin() const noexcept
{
    return alleles_.cbegin();
}

typename Genotype<Allele>::Iterator Genotype<Allele>::cend() const noexcept
{
    return alleles_.cend();
}

bool Genotype<Allele>::is_homozygous() const
{
    return std::adjacent_find(std::cbegin(alleles_), std::cend(alleles_), std::not_equal_to<Allele>()) == std::cend(alleles_);
}

unsigned Genotype<Allele>::zygosity() const
{
    return static_cast<unsigned>(get_unique().size());
}

bool Genotype<Allele>::contains(const Allele& allele) const
{
    return std::find(std::cbegin(alleles_), std::cend(alleles_), allele) != std::cend(alleles_);
}

unsigned Genotype<Allele>::count(const Allele& element) const
{
    return static_cast<unsigned>(std::count(std::cbegin(alleles_), std::cend(alleles_), element));
}

unsigned Genotype<Allele>::ploidy() const noexcept
{
    return static_cast<unsigned>(alleles_.size());
}

std::vector<Allele> Genotype<Allele>::get_unique() const
{
    std::vector<Allele> result {};
    result.reserve(ploidy());
    
    std::unique_copy(std::cbegin(alleles_), std::cend(alleles_), std::back_inserter(result));
    
    return result;
}

// non-member methods

bool is_homozygous_reference(const Genotype<Haplotype>& genotype, const Allele& reference)
{
    return splice<Allele>(genotype, get_region(reference)).count(reference) == genotype.ploidy();
}

unsigned num_genotypes(unsigned num_elements, unsigned ploidy)
{
    return static_cast<unsigned>(boost::math::binomial_coefficient<double>(num_elements + ploidy - 1, num_elements - 1));
}

void print_alleles(const Genotype<Haplotype>& genotype)
{
    for (unsigned i {}; i < genotype.ploidy() - 1; ++i) {
        print_alleles(genotype[i]);
        std::cout << std::endl;
    }
    print_alleles(genotype[genotype.ploidy() - 1]);
}
