//
//  cancer_genotype.cpp
//  Octopus
//
//  Created by Daniel Cooke on 09/04/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#include "cancer_genotype.hpp"

#include <cassert>
#include <iostream>

namespace Octopus
{
std::size_t estimate_num_cancer_genotypes(const std::size_t num_haplotypes, const unsigned ploidy)
{
    return (num_genotypes(static_cast<unsigned>(num_haplotypes), ploidy) - 1) * num_haplotypes;
}

namespace
{
    auto make_all_shared(const std::vector<Haplotype>& elements)
    {
        std::vector<std::shared_ptr<Haplotype>> result(elements.size());
        std::transform(std::cbegin(elements), std::cend(elements), std::begin(result),
                       [] (const auto& element) { return std::make_shared<Haplotype>(element); });
        return result;
    }
} // namespace

std::vector<CancerGenotype<Haplotype>>
generate_all_cancer_genotypes(const std::vector<Haplotype>& haplotypes, const unsigned ploidy)
{
    assert(!haplotypes.empty());
    
    const auto haplotypes_ptrs = make_all_shared(haplotypes);
    
    const auto germline_genotypes = generate_all_genotypes(haplotypes_ptrs, ploidy);
    
    std::vector<CancerGenotype<Haplotype>> result {};
    
    if (haplotypes.size() > ploidy) {
        result.reserve(estimate_num_cancer_genotypes(haplotypes.size(), ploidy));
        for (const auto& germline_genotype : germline_genotypes) {
            for (const auto& ptr : haplotypes_ptrs) {
                if (!contains(germline_genotype, *ptr)) {
                    result.emplace_back(germline_genotype, ptr);
                }
            }
        }
    } else if (ploidy == 1) {
        result.reserve(1);
        result.emplace_back(germline_genotypes.front(), haplotypes_ptrs.front());
    } else if (ploidy == 2) {
        if (haplotypes.size() == 1) {
            result.reserve(1);
            result.emplace_back(germline_genotypes.front(), haplotypes_ptrs.front());
        } else {
            assert(!germline_genotypes[1].is_homozygous());
            assert(germline_genotypes[0].contains(haplotypes[0]));
            result.reserve(3);
            result.emplace_back(germline_genotypes[0], haplotypes_ptrs[1]);
            result.emplace_back(germline_genotypes[1], haplotypes_ptrs[0]);
            result.emplace_back(germline_genotypes[2], haplotypes_ptrs[0]);
        }
    }
    else {
        result.reserve(germline_genotypes.size());
        for (const auto& germline_genotype : germline_genotypes) {
            for (const auto& ptr : haplotypes_ptrs) {
                if (germline_genotype.is_homozygous() || !contains(germline_genotype, *ptr)) {
                    result.emplace_back(germline_genotype, ptr);
                }
            }
        }
    }
    
    result.shrink_to_fit();
    
    return result;
}

namespace debug
{
    void print_alleles(const CancerGenotype<Haplotype>& genotype)
    {
        print_alleles(std::cout, genotype);
    }
    
    void print_variant_alleles(const CancerGenotype<Haplotype>& genotype)
    {
        print_variant_alleles(std::cout, genotype);
    }
} // namespace debug
} // namespace Octopus