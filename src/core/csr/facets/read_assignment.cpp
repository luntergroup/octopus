// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "read_assignment.hpp"

#include <utility>
#include <iterator>
#include <algorithm>
#include <limits>
#include <iostream>
#include <stdexcept>

#include <utils/maths.hpp>
#include <utils/append.hpp>
#include <core/models/haplotype_likelihood_model.hpp>
#include <core/models/utils/kmer_mapper.hpp>

namespace octopus {

using HaplotypeLikelihoods = std::vector<std::vector<double>>;

auto max_posterior_haplotypes(const Genotype<Haplotype>& genotype, const unsigned read,
                              const HaplotypeLikelihoods& likelihoods)
{
    std::vector<unsigned> result {};
    result.reserve(genotype.ploidy());
    
    auto max_likelihood = std::numeric_limits<double>::lowest();
    
    for (unsigned k {0}; k < genotype.ploidy(); ++k) {
        const auto curr = likelihoods[k][read];
        
        if (maths::almost_equal(curr, max_likelihood)) {
            result.push_back(k);
        } else if (curr > max_likelihood) {
            result.assign({k});
            max_likelihood = curr;
        }
    }
    
    result.erase(std::unique(std::begin(result), std::end(result)), std::end(result));
    
    return result;
}

auto calculate_support(const Genotype<Haplotype>& genotype,
                       const std::vector<AlignedRead>& reads,
                       const HaplotypeLikelihoods& likelihoods)
{
    HaplotypeSupportMap result {};
    
    for (unsigned i {0}; i < reads.size(); ++i) {
        const auto top = max_posterior_haplotypes(genotype, i, likelihoods);
        
        if (top.size() == 1) {
            result[genotype[top.front()]].push_back(reads[i]);
        }
    }
    
    return result;
}

auto expand(const Genotype<Haplotype>& genotype, unsigned n)
{
    Genotype<Haplotype> result {genotype.ploidy()};
    
    for (const auto& haplotype : genotype) {
        result.emplace(expand(haplotype, n));
    }
    
    return result;
}

auto calculate_likelihood(const Haplotype& haplotype,
                          const AlignedRead& read,
                          const HaplotypeLikelihoodModel& model)
{
    static constexpr unsigned char K {6};
    auto mapping_positions = map_query_to_target<K>(read.sequence(), haplotype.sequence());
    return model.ln_probability(read, mapping_positions);
}

auto calculate_likelihoods(const Haplotype& haplotype,
                           const std::vector<AlignedRead>& reads,
                           HaplotypeLikelihoodModel& model)
{
    model.reset(haplotype);
    HaplotypeLikelihoods::value_type result(reads.size());
    std::transform(std::cbegin(reads), std::cend(reads), std::begin(result),
                   [&] (const auto& read) {
                       return calculate_likelihood(haplotype, read, model);
                   });
    return result;
}

auto calculate_likelihoods(const Genotype<Haplotype>& genotype,
                           const std::vector<AlignedRead>& reads,
                           HaplotypeLikelihoodModel& model)
{
    const auto& genotype_region = mapped_region(genotype);
    
    const auto reads_region = encompassing_region(reads);
    
    unsigned min_lhs_expansion {0};
    
    if (begins_before(reads_region, genotype_region)) {
        min_lhs_expansion += begin_distance(reads_region, genotype_region);
    }
    
    unsigned min_rhs_expansion {0};
    
    if (ends_before(genotype_region, reads_region)) {
        min_rhs_expansion += end_distance(genotype_region, reads_region);
    }
    
    const auto min_expansion = std::max({min_lhs_expansion, min_rhs_expansion, 20u});
    
    std::map<Haplotype, Haplotype> expanded_haplotypes {};
    
    for (const auto& haplotype : genotype) {
        expanded_haplotypes.emplace(haplotype, expand(haplotype, min_expansion));
    }
    
    const auto expanded_genotype = expand(genotype, min_expansion);
    
    const auto haplotypes = expanded_genotype.copy_unique();
    
    HaplotypeLikelihoods result(genotype.ploidy());
    
    std::transform(std::cbegin(genotype), std::cend(genotype), std::begin(result),
                   [&] (const auto& haplotype) {
                       return calculate_likelihoods(expanded_haplotypes.at(haplotype), reads, model);
                   });
    
    return result;
}

HaplotypeSupportMap compute_haplotype_support(const Genotype<Haplotype>& genotype,
                                              const std::vector<AlignedRead>& reads)
{
    return compute_haplotype_support(genotype, reads, HaplotypeLikelihoodModel {});
}

HaplotypeSupportMap compute_haplotype_support(const Genotype<Haplotype>& genotype,
                                              const std::vector<AlignedRead>& reads,
                                              HaplotypeLikelihoodModel model)
{
    const auto likelihoods = calculate_likelihoods(genotype, reads, model);
    return calculate_support(genotype, reads, likelihoods);
}

AlleleSupportMap compute_allele_support(const std::vector<Allele>& alleles,
                                        const HaplotypeSupportMap& haplotype_support)
{
    AlleleSupportMap result {};
    result.reserve(alleles.size());
    
    for (const auto& allele : alleles) {
        for (const auto& p : haplotype_support) {
            if (p.first.contains(allele)) {
                utils::append(p.second, result[allele]);
            }
        }
    }
    
    return result;
    
}

} // namespace octopus
