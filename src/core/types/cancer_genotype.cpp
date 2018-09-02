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

std::vector<CancerGenotype<Haplotype>>
generate_all_cancer_genotypes(const std::vector<Genotype<Haplotype>>& germline_genotypes,
                              const std::vector<Haplotype>& somatic_haplotypes,
                              const unsigned somatic_ploidy, const bool allow_shared)
{
    const auto somatic_genotypes = generate_all_max_zygosity_genotypes(somatic_haplotypes, somatic_ploidy);
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
    const auto somatic_genotypes = generate_all_max_zygosity_genotypes(somatic_haplotypes, somatic_ploidy,
                                                                       somatic_genotype_indices);
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

std::vector<CancerGenotype<Haplotype>>
extend_somatic_genotypes(const std::vector<CancerGenotype<Haplotype>>& old_genotypes,
                         const std::vector<Haplotype>& somatic_haplotypes,
                         const bool allow_shared)
{
    std::vector<std::shared_ptr<Haplotype>> temp_pointers(somatic_haplotypes.size());
    std::transform(std::cbegin(somatic_haplotypes), std::cend(somatic_haplotypes), std::begin(temp_pointers),
                   [] (const auto& haplotype) { return std::make_shared<Haplotype>(haplotype); });
    std::vector<CancerGenotype<Haplotype>> result {};
    result.reserve(old_genotypes.size() * somatic_haplotypes.size());
    for (const auto& old_genotype : old_genotypes) {
        for (const auto& haplotype_ptr : temp_pointers) {
            if (allow_shared || !old_genotype.germline().contains(*haplotype_ptr)) {
                auto new_somatic_genotype = old_genotype.somatic();
                new_somatic_genotype.emplace(haplotype_ptr);
                if (is_max_zygosity(new_somatic_genotype)) {
                    result.emplace_back(old_genotype.germline(), std::move(new_somatic_genotype));
                }
            }
        }
    }
    return result;
}

std::vector<CancerGenotype<Haplotype>>
extend_somatic_genotypes(const std::vector<CancerGenotype<Haplotype>>& old_genotypes,
                         const std::vector<CancerGenotypeIndex>& old_genotype_indices,
                         const std::vector<Haplotype>& somatic_haplotypes,
                         std::vector<CancerGenotypeIndex>& new_genotype_indices,
                         const bool allow_shared)
{
    assert(old_genotypes.size() == old_genotype_indices.size());
    std::vector<std::shared_ptr<Haplotype>> temp_pointers(somatic_haplotypes.size());
    std::transform(std::cbegin(somatic_haplotypes), std::cend(somatic_haplotypes), std::begin(temp_pointers),
                   [] (const auto& haplotype) { return std::make_shared<Haplotype>(haplotype); });
    std::vector<CancerGenotype<Haplotype>> result {};
    const auto max_new_genotypes = old_genotypes.size() * somatic_haplotypes.size();
    result.reserve(max_new_genotypes);
    new_genotype_indices.clear();
    new_genotype_indices.reserve(max_new_genotypes);
    for (std::size_t g {0}; g < old_genotypes.size(); ++g) {
        const auto& old_genotype = old_genotypes[g];
        for (std::size_t h {0}; h < somatic_haplotypes.size(); ++h) {
            const auto& haplotype_ptr = temp_pointers[h];
            if (allow_shared || !old_genotype.germline().contains(*haplotype_ptr)) {
                auto new_somatic_genotype = old_genotype.somatic();
                new_somatic_genotype.emplace(haplotype_ptr);
                if (is_max_zygosity(new_somatic_genotype)) {
                    result.emplace_back(old_genotype.germline(), std::move(new_somatic_genotype));
                    auto new_genotype_index = old_genotype_indices[g];
                    new_genotype_index.somatic.push_back(h);
                    new_genotype_indices.push_back(std::move(new_genotype_index));
                }
            }
        }
    }
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
