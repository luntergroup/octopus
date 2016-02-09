//
//  genotype.cpp
//  Octopus
//
//  Created by Daniel Cooke on 02/10/2015.
//  Copyright © 2015 Oxford University. All rights reserved.
//

#include "genotype.hpp"

#include <iostream>

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

const GenomicRegion& Genotype<Allele>::get_region() const noexcept
{
    return alleles_.front().get_region();
}

const Allele& Genotype<Allele>::at(unsigned n) const
{
    return alleles_.at(n);
}

const Allele& Genotype<Allele>::operator[](unsigned n) const
{
    return alleles_[n];
}

unsigned Genotype<Allele>::ploidy() const noexcept
{
    return static_cast<unsigned>(alleles_.size());
}

bool Genotype<Allele>::is_homozygous() const
{
    return std::adjacent_find(std::cbegin(alleles_), std::cend(alleles_),
                              std::not_equal_to<Allele>()) == std::cend(alleles_);
}

unsigned Genotype<Allele>::zygosity() const
{
    if (ploidy() == 1 || is_homozygous()) {
        return 1;
    } else if (ploidy() == 2) {
        return 2;
    }
    return static_cast<unsigned>(copy_unique().size());
}

bool Genotype<Allele>::contains(const Allele& allele) const
{
    return std::find(std::cbegin(alleles_), std::cend(alleles_), allele) != std::cend(alleles_);
}

unsigned Genotype<Allele>::count(const Allele& element) const
{
    return static_cast<unsigned>(std::count(std::cbegin(alleles_), std::cend(alleles_), element));
}

std::vector<Allele> Genotype<Allele>::copy_unique() const
{
    auto result = alleles_;
    
    std::sort(std::begin(result), std::end(result));
    
    result.erase(std::unique(std::begin(result), std::end(result)), std::end(result));
    
    return result;
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

// non-member methods

bool contains(const Genotype<Haplotype>& genotype, const Allele& allele)
{
    return std::any_of(std::cbegin(genotype), std::cend(genotype),
                       [&allele] (const auto& haplotype) { return haplotype.contains(allele); });
}

bool contains_exact(const Genotype<Haplotype>& genotype, const Allele& allele)
{
    return std::any_of(std::cbegin(genotype), std::cend(genotype),
                       [&allele] (const auto& haplotype) { return haplotype.contains_exact(allele); });
}

bool is_homozygous(const Genotype<Haplotype>& genotype, const Allele& allele)
{
    return splice<Allele>(genotype, get_region(allele)).count(allele) == genotype.ploidy();
}

size_t num_genotypes(const unsigned num_elements, const unsigned ploidy)
{
    return static_cast<size_t>(boost::math::binomial_coefficient<double>(num_elements + ploidy - 1,
                                                                         num_elements - 1));
}

void print_alleles(const Genotype<Haplotype>& genotype)
{
    if (genotype.ploidy() == 0) {
        std::cout << "[]";
    }
    const auto element_counts = make_element_count_map(genotype);
    std::vector<std::pair<Haplotype, unsigned>> p {element_counts.begin(), element_counts.end()};
    std::cout << "[";
    for (unsigned i {}; i < p.size() - 1; ++i) {
        print_alleles(p[i].first);
        std::cout << "(" << p[i].second << "),";
    }
    print_alleles(p.back().first);
    std::cout << "(" << p.back().second << ")]";
}

void print_variant_alleles(const Genotype<Haplotype>& genotype)
{
    if (genotype.ploidy() == 0) {
        std::cout << "[]";
    }
    const auto element_counts = make_element_count_map(genotype);
    std::vector<std::pair<Haplotype, unsigned>> p {element_counts.begin(), element_counts.end()};
    std::cout << "[";
    for (unsigned i {}; i < p.size() - 1; ++i) {
        print_variant_alleles(p[i].first);
        std::cout << "(" << p[i].second << "),";
    }
    print_variant_alleles(p.back().first);
    std::cout << "(" << p.back().second << ")]";
}
