// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "cancer_genotype.hpp"

#include <cassert>
#include <iostream>

namespace octopus {

bool contains(const CancerGenotype<Haplotype>& genotype, const Allele& allele)
{
    return contains(genotype.germline(), allele) || contains(genotype.somatic(), allele);
}

bool includes(const CancerGenotype<Haplotype>& genotype, const Allele& allele)
{
    return includes(genotype.germline(), allele) || includes(genotype.somatic(), allele);
}

namespace {

auto make_all_shared(const std::vector<Haplotype>& elements)
{
    std::vector<std::shared_ptr<Haplotype>> result(elements.size());
    std::transform(std::cbegin(elements), std::cend(elements), std::begin(result),
                   [] (const auto& element) { return std::make_shared<Haplotype>(element); });
    return result;
}

bool is_proper_somatic_genotype(const Genotype<Haplotype>& genotype)
{
    return genotype.zygosity() == genotype.ploidy();
}

std::vector<Genotype<Haplotype>> generate_somatic_genotypes(const std::vector<Haplotype>& haplotypes,
                                                            const unsigned ploidy)
{
    auto result = generate_all_genotypes(haplotypes, ploidy);
    if (ploidy > 1) {
        auto itr = std::remove_if(std::begin(result), std::end(result), is_proper_somatic_genotype);
        result.erase(itr, std::end(result));
    }
    return result;
}

std::vector<Genotype<Haplotype>> generate_somatic_genotypes(const std::vector<Haplotype>& haplotypes,
                                                            const unsigned ploidy,
                                                            std::vector<GenotypeIndex>& germline_genotype_indices)
{
    std::vector<GenotypeIndex> all_genotype_indices {};
    auto all_genotypes = generate_all_genotypes(haplotypes, ploidy, all_genotype_indices);
    if (ploidy > 1) {
        std::vector<Genotype<Haplotype>> result {};
        result.reserve(all_genotypes.size());
        for (std::size_t i {0}; i < all_genotypes.size(); ++i) {
            if (is_proper_somatic_genotype(all_genotypes[i])) {
                result.push_back(std::move(all_genotypes[i]));
                germline_genotype_indices.push_back(std::move(all_genotype_indices[i]));
            }
        }
        return result;
    } else {
        germline_genotype_indices = std::move(all_genotype_indices);
        return all_genotypes;
    }
}

} // namespace

std::vector<CancerGenotype<Haplotype>>
generate_all_cancer_genotypes(const std::vector<Genotype<Haplotype>>& germline_genotypes,
                              const std::vector<Haplotype>& somatic_haplotypes,
                              const unsigned somatic_ploidy, const bool allow_shared)
{
    const auto somatic_genotypes = generate_somatic_genotypes(somatic_haplotypes, somatic_ploidy);
    std::vector<CancerGenotype<Haplotype>> result {};
    result.reserve(germline_genotypes.size() * somatic_genotypes.size());
    for (const auto& germline : germline_genotypes) {
        for (const auto& somatic : somatic_genotypes) {
            if (allow_shared || !have_shared(germline, somatic)) {
                result.emplace_back(germline, somatic);
            }
        }
    }
    return result;
}

std::vector<CancerGenotype<Haplotype>>
generate_all_cancer_genotypes(const std::vector<Genotype<Haplotype>>& germline_genotypes,
                              const std::vector<GenotypeIndex>& germline_genotype_indices,
                              const std::vector<Haplotype>& somatic_haplotypes,
                              std::vector<CancerGenotypeIndex>& cancer_genotype_indices,
                              const unsigned somatic_ploidy, const bool allow_shared)
{
    assert(germline_genotypes.size() == germline_genotype_indices.size());
    std::vector<GenotypeIndex> somatic_genotype_indices {};
    const auto somatic_genotypes = generate_somatic_genotypes(somatic_haplotypes, somatic_ploidy, somatic_genotype_indices);
    std::vector<CancerGenotype<Haplotype>> result {};
    const auto max_cancer_genotypes = germline_genotypes.size() * somatic_genotypes.size();
    result.reserve(max_cancer_genotypes);
    cancer_genotype_indices.reserve(max_cancer_genotypes);
    for (std::size_t g {0}; g < germline_genotypes.size(); ++g) {
        for (std::size_t h {0}; h < somatic_genotypes.size(); ++h) {
            if (allow_shared || !have_shared(germline_genotypes[g], somatic_genotypes[h])) {
                result.emplace_back(germline_genotypes[g], somatic_genotypes[h]);
                cancer_genotype_indices.push_back({germline_genotype_indices[g], somatic_genotype_indices[h]});
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
