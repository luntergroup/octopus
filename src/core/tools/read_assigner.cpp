// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "read_assigner.hpp"

#include <utility>
#include <iterator>
#include <algorithm>
#include <limits>
#include <stdexcept>
#include <cassert>

#include "utils/maths.hpp"
#include "core/models/haplotype_likelihood_model.hpp"
#include "utils/kmer_mapper.hpp"

namespace octopus {

using HaplotypeLikelihoods = std::vector<std::vector<double>>;

auto max_posterior_haplotypes(const std::vector<Haplotype>& haplotypes, const unsigned read,
                              const HaplotypeLikelihoods& likelihoods)
{
    std::vector<unsigned> result {};
    result.reserve(haplotypes.size());
    auto max_likelihood = std::numeric_limits<double>::lowest();
    for (unsigned k {0}; k < haplotypes.size(); ++k) {
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

auto calculate_support(const std::vector<Haplotype>& haplotypes,
                       const std::vector<AlignedRead>& reads,
                       const HaplotypeLikelihoods& likelihoods)
{
    HaplotypeSupportMap result {};
    for (unsigned i {0}; i < reads.size(); ++i) {
        const auto top = max_posterior_haplotypes(haplotypes, i, likelihoods);
        if (top.size() == 1) {
            result[haplotypes[top.front()]].push_back(reads[i]);
        }
    }
    return result;
}

auto max_deletion_size(const std::vector<Haplotype>& haplotypes)
{
    unsigned result {0};
    for (const auto& haplotype : haplotypes) {
        if (region_size(haplotype) > sequence_size(haplotype)) {
            auto diff = static_cast<unsigned>(region_size(haplotype) - sequence_size(haplotype));
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

auto calculate_likelihoods(const std::vector<Haplotype>& haplotypes,
                           const std::vector<AlignedRead>& reads,
                           HaplotypeLikelihoodModel& model)
{
    assert(!haplotypes.empty());
    const auto& genotype_region = mapped_region(haplotypes.front());
    const auto reads_region = encompassing_region(reads);
    const auto min_flank_pad = HaplotypeLikelihoodModel::pad_requirement();
    unsigned min_lhs_expansion {min_flank_pad}, min_rhs_expansion {min_flank_pad};
    if (begins_before(reads_region, genotype_region)) {
        min_lhs_expansion += begin_distance(reads_region, genotype_region);
    }
    if (ends_before(genotype_region, reads_region)) {
        min_rhs_expansion += end_distance(genotype_region, reads_region);
    }
    const auto min_expansion = std::max({min_lhs_expansion, min_rhs_expansion, 20u}) + max_deletion_size(haplotypes);
    const auto read_hashes = compute_read_hashes(reads);
    static constexpr unsigned char mapperKmerSize {6};
    auto haplotype_hashes = init_kmer_hash_table<mapperKmerSize>();
    HaplotypeLikelihoods result {};
    result.reserve(haplotypes.size());
    for (const auto& haplotype : haplotypes) {
        const auto expanded_haplotype = expand(haplotype, min_expansion);
        populate_kmer_hash_table<mapperKmerSize>(expanded_haplotype.sequence(), haplotype_hashes);
        auto haplotype_mapping_counts = init_mapping_counts(haplotype_hashes);
        model.reset(expanded_haplotype);
        std::vector<double> likelihoods(reads.size());
        std::transform(std::cbegin(reads), std::cend(reads), std::cbegin(read_hashes), std::begin(likelihoods),
                       [&] (const auto& read, const auto& read_hash) {
                           auto mapping_positions = map_query_to_target(read_hash, haplotype_hashes, haplotype_mapping_counts);
                           reset_mapping_counts(haplotype_mapping_counts);
                           return model.evaluate(read, mapping_positions);
                       });
        clear_kmer_hash_table(haplotype_hashes);
        result.push_back(std::move(likelihoods));
    }
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
        const auto unique_haplotypes = genotype.copy_unique();
        assert(unique_haplotypes.size() > 1);
        const auto likelihoods = calculate_likelihoods(unique_haplotypes, reads, model);
        return calculate_support(unique_haplotypes, reads, likelihoods);
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
