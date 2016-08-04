// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "cancer_genotype.hpp"

#include <cassert>
#include <iostream>

namespace octopus {

bool contains(const CancerGenotype<Haplotype>& genotype, const Allele& allele)
{
    return contains(genotype.germline_genotype(), allele) || genotype.somatic_element().contains(allele);
}

bool includes(const CancerGenotype<Haplotype>& genotype, const Allele& allele)
{
    return includes(genotype.germline_genotype(), allele) || genotype.somatic_element().includes(allele);
}
    
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

std::pair<std::vector<CancerGenotype<Haplotype>>, std::vector<Genotype<Haplotype>>>
generate_all_cancer_genotypes(const std::vector<Haplotype>& haplotypes, const unsigned ploidy)
{
    assert(!haplotypes.empty());
    
    const auto haplotypes_ptrs = make_all_shared(haplotypes);
    
    const auto germline_genotypes = generate_all_genotypes(haplotypes_ptrs, ploidy);
    
    std::vector<CancerGenotype<Haplotype>> cancer_genotypes {};
    cancer_genotypes.reserve(estimate_num_cancer_genotypes(haplotypes.size(), ploidy));
    
    for (const auto& germline_genotype : germline_genotypes) {
        for (const auto& ptr : haplotypes_ptrs) {
            if (!contains(germline_genotype, *ptr)) {
                cancer_genotypes.emplace_back(germline_genotype, ptr);
            }
        }
    }
    
    cancer_genotypes.shrink_to_fit();
    
    return std::make_pair(std::move(cancer_genotypes), std::move(germline_genotypes));
}

namespace debug {
    void print_alleles(const CancerGenotype<Haplotype>& genotype)
    {
        print_alleles(std::cout, genotype);
    }
    
    void print_variant_alleles(const CancerGenotype<Haplotype>& genotype)
    {
        print_variant_alleles(std::cout, genotype);
    }
} // namespace debug
} // namespace octopus
