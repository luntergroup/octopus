// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "read_assigner.hpp"

#include <utility>
#include <iterator>
#include <algorithm>
#include <limits>
#include <iostream>
#include <stdexcept>

#include "utils/maths.hpp"
#include "core/models/haplotype_likelihood_model.hpp"
#include "utils/kmer_mapper.hpp"

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

auto max_deletion_size(const Genotype<Haplotype>& genotype)
{
    unsigned result {0};
    for (const auto& haplotype : genotype) {
        if (region_size(genotype) > sequence_size(haplotype)) {
            auto diff = static_cast<unsigned>(region_size(genotype) - sequence_size(haplotype));
            result = std::max(result, diff);
        }
    }
    return result;
}

auto compute_read_hashes(const std::vector<AlignedRead>& reads)
{
    static constexpr unsigned char mapperKmerSize {6};
    std::vector<KmerPerfectHashes> result {};
    result.reserve(reads.size());
    std::transform(std::cbegin(reads), std::cend(reads), std::back_inserter(result),
                   [=] (const AlignedRead& read) { return compute_kmer_hashes<mapperKmerSize>(read.sequence()); });
    return result;
}

auto calculate_likelihoods(const Genotype<Haplotype>& genotype,
                           const std::vector<AlignedRead>& reads,
                           HaplotypeLikelihoodModel& model)
{
    const auto& genotype_region = mapped_region(genotype);
    const auto reads_region = encompassing_region(reads);
    const auto min_flank_pad = HaplotypeLikelihoodModel::pad_requirement();
    unsigned min_lhs_expansion {min_flank_pad}, min_rhs_expansion {min_flank_pad};
    if (begins_before(reads_region, genotype_region)) {
        min_lhs_expansion += begin_distance(reads_region, genotype_region);
    }
    if (ends_before(genotype_region, reads_region)) {
        min_rhs_expansion += end_distance(genotype_region, reads_region);
    }
    const auto min_expansion = std::max({min_lhs_expansion, min_rhs_expansion, 20u}) + max_deletion_size(genotype);
    std::map<Haplotype, Haplotype> expanded_haplotypes {};
    for (const auto& haplotype : genotype) {
        expanded_haplotypes.emplace(haplotype, expand(haplotype, min_expansion));
    }
    const auto read_hashes = compute_read_hashes(reads);
    static constexpr unsigned char mapperKmerSize {6};
    auto haplotype_hashes = init_kmer_hash_table<mapperKmerSize>();
    std::map<Haplotype, std::vector<double>> likelihoods {};
    for (const auto& p : expanded_haplotypes) {
        populate_kmer_hash_table<mapperKmerSize>(p.second.sequence(), haplotype_hashes);
        auto haplotype_mapping_counts = init_mapping_counts(haplotype_hashes);
        model.reset(p.second);
        likelihoods[p.first].resize(reads.size());
        std::transform(std::cbegin(reads), std::cend(reads), std::cbegin(read_hashes), std::begin(likelihoods[p.first]),
                       [&] (const auto& read, const auto& read_hash) {
                           auto mapping_positions = map_query_to_target(read_hash, haplotype_hashes, haplotype_mapping_counts);
                           reset_mapping_counts(haplotype_mapping_counts);
                           return model.evaluate(read, mapping_positions);
                       });
        clear_kmer_hash_table(haplotype_hashes);
    }
    HaplotypeLikelihoods result(genotype.ploidy());
    std::transform(std::cbegin(genotype), std::cend(genotype), std::begin(result),
                   [&] (const auto& haplotype) { return likelihoods[haplotype]; });
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
    if (!genotype.is_homozygous() && !reads.empty()) {
        const auto likelihoods = calculate_likelihoods(genotype, reads, model);
        return calculate_support(genotype, reads, likelihoods);
    } else {
        return {};
    }
}

AlleleSupportMap compute_allele_support(const std::vector<Allele>& alleles,
                                        const HaplotypeSupportMap& haplotype_support)
{
    AlleleSupportMap result {};
    result.reserve(alleles.size());
    for (const auto& allele : alleles) {
        ReadRefSupportSet allele_support {};
        for (const auto& p : haplotype_support) {
            if (p.first.contains(allele)) {
                allele_support.insert(std::cend(allele_support), std::cbegin(p.second), std::cend(p.second));
            }
        }
        result.emplace(allele, std::move(allele_support));
    }
    return result;
}

} // namespace octopus
