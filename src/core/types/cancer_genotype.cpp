// Copyright (c) 2015-2018 Daniel Cooke
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

namespace {

auto make_all_shared(const std::vector<Haplotype>& elements)
{
    std::vector<std::shared_ptr<Haplotype>> result(elements.size());
    std::transform(std::cbegin(elements), std::cend(elements), std::begin(result),
                   [] (const auto& element) { return std::make_shared<Haplotype>(element); });
    return result;
}

} // namespace

std::vector<CancerGenotype<Haplotype>>
generate_all_cancer_genotypes(const std::vector<Genotype<Haplotype>>& germline_genotypes,
                              const std::vector<Haplotype>& somatic_haplotypes)
{
    const auto haplotype_ptrs = make_all_shared(somatic_haplotypes);
    std::vector<CancerGenotype<Haplotype>> result {};
    result.reserve(germline_genotypes.size() * somatic_haplotypes.size());
    for (const auto& genotype : germline_genotypes) {
        for (const auto& haplotype : haplotype_ptrs) {
            if (!contains(genotype, *haplotype)) {
                result.emplace_back(genotype, haplotype);
            }
        }
    }
    return result;
}

std::vector<CancerGenotype<Haplotype>>
generate_all_cancer_genotypes(const std::vector<Genotype<Haplotype>>& germline_genotypes,
                              const std::vector<std::vector<unsigned>>& germline_genotype_indices,
                              const std::vector<Haplotype>& somatic_haplotypes,
                              std::vector<std::pair<std::vector<unsigned>, unsigned>>& cancer_genotype_indices)
{
    assert(germline_genotypes.size() == germline_genotype_indices.size());
    const auto haplotype_ptrs = make_all_shared(somatic_haplotypes);
    std::vector<CancerGenotype<Haplotype>> result {};
    const auto max_cancer_genotypes = germline_genotypes.size() * somatic_haplotypes.size();
    result.reserve(max_cancer_genotypes);
    cancer_genotype_indices.reserve(max_cancer_genotypes);
    for (std::size_t g {0}; g < germline_genotypes.size(); ++g) {
        for (unsigned h {0}; h < somatic_haplotypes.size(); ++h) {
            if (!contains(germline_genotypes[g], *haplotype_ptrs[h])) {
                result.emplace_back(germline_genotypes[g], haplotype_ptrs[h]);
                cancer_genotype_indices.emplace_back(germline_genotype_indices[g], h);
            }
        }
    }
    cancer_genotype_indices.shrink_to_fit();
    return result;
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
