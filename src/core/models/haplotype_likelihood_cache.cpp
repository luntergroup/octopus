// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "haplotype_likelihood_cache.hpp"

#include <utility>
#include <cassert>

#include <iostream> // DEBUG
#include <iomanip>  // DEBUG

namespace octopus {

// public methods

HaplotypeLikelihoodCache::HaplotypeLikelihoodCache(const unsigned max_haplotypes,
                                                   const std::vector<SampleName>& samples)
: cache_ {max_haplotypes}
, sample_indices_ {samples.size()}
{
    mapping_positions_.resize(maxMappingPositions);
}

HaplotypeLikelihoodCache::HaplotypeLikelihoodCache(HaplotypeLikelihoodModel likelihood_model,
                                                   unsigned max_haplotypes,
                                                   const std::vector<SampleName>& samples)
: likelihood_model_ {std::move(likelihood_model)}
, cache_ {max_haplotypes}
, sample_indices_ {samples.size()}
{
    mapping_positions_.resize(maxMappingPositions);
}

HaplotypeLikelihoodCache::ReadPacket::ReadPacket(Iterator first, Iterator last)
: first {first}
, last {last}
, num_reads {static_cast<std::size_t>(std::distance(first, last))}
{}

void HaplotypeLikelihoodCache::populate(const ReadMap& reads,
                                        const std::vector<Haplotype>& haplotypes,
                                        boost::optional<FlankState> flank_state)
{
    // This code is not very pretty because it is a bottleneck for the entire application.
    // We want to try a minimise memory allocations for the mapping.
    cache_.clear();
    if (cache_.bucket_count() < haplotypes.size()) {
        cache_.rehash(haplotypes.size());
    }
    set_read_iterators_and_sample_indices(reads);
    assert(reads.size() == read_iterators_.size());
    const auto num_samples = reads.size();
    // Precompute all read hashes so we don't have to recompute for each haplotype
    std::vector<std::vector<KmerPerfectHashes>> read_hashes {};
    read_hashes.reserve(num_samples);
    for (const auto& t : read_iterators_) {
        std::vector<KmerPerfectHashes> sample_read_hashes {};
        sample_read_hashes.reserve(t.num_reads);
        std::transform(t.first, t.last, std::back_inserter(sample_read_hashes),
                       [] (const AlignedRead& read) { return compute_kmer_hashes<mapperKmerSize>(read.sequence()); });
        read_hashes.emplace_back(std::move(sample_read_hashes));
    }
    auto haplotype_hashes = init_kmer_hash_table<mapperKmerSize>();
    const auto first_mapping_position = std::begin(mapping_positions_);
    for (const auto& haplotype : haplotypes) {
        populate_kmer_hash_table<mapperKmerSize>(haplotype.sequence(), haplotype_hashes);
        auto haplotype_mapping_counts = init_mapping_counts(haplotype_hashes);
        auto itr = std::begin(cache_.emplace(std::piecewise_construct,
                                             std::forward_as_tuple(haplotype),
                                             std::forward_as_tuple(num_samples)).first->second);
        likelihood_model_.reset(haplotype, flank_state);
        auto read_hash_itr = std::cbegin(read_hashes);
        for (const auto& t : read_iterators_) { // for each sample
            *itr = std::vector<double>(t.num_reads);
            std::transform(t.first, t.last, std::cbegin(*read_hash_itr), std::begin(*itr),
                           [&] (const AlignedRead& read, const auto& read_hashes) {
                               const auto last_mapping_position = map_query_to_target(read_hashes, haplotype_hashes,
                                                                                      haplotype_mapping_counts,
                                                                                      first_mapping_position,
                                                                                      maxMappingPositions);
                               reset_mapping_counts(haplotype_mapping_counts);
                               return likelihood_model_.evaluate(read, first_mapping_position, last_mapping_position);
                           });
            ++read_hash_itr;
            ++itr;
        }
        clear_kmer_hash_table(haplotype_hashes);
    }
    likelihood_model_.clear();
    read_iterators_.clear();
}

std::size_t HaplotypeLikelihoodCache::num_likelihoods(const SampleName& sample) const
{
    return std::cbegin(cache_)->second.at(sample_indices_.at(sample)).size();
}

const HaplotypeLikelihoodCache::LikelihoodVector&
HaplotypeLikelihoodCache::operator()(const SampleName& sample, const Haplotype& haplotype) const
{
    return cache_.at(haplotype)[sample_indices_.at(sample)];
}

const HaplotypeLikelihoodCache::LikelihoodVector&
HaplotypeLikelihoodCache::operator[](const Haplotype& haplotype) const
{
    return cache_.at(haplotype)[*primed_sample_];
}

HaplotypeLikelihoodCache::SampleLikelihoodMap
HaplotypeLikelihoodCache::extract_sample(const SampleName& sample) const
{
    const auto sample_index = sample_indices_.at(sample);
    SampleLikelihoodMap result {cache_.size()};
    for (const auto& p : cache_) {
        result.emplace(p.first, p.second[sample_index]);
    }
    return result;
}

bool HaplotypeLikelihoodCache::contains(const Haplotype& haplotype) const noexcept
{
    return cache_.count(haplotype) == 1;
}

bool HaplotypeLikelihoodCache::is_empty() const noexcept
{
    return cache_.empty();
}

void HaplotypeLikelihoodCache::clear() noexcept
{
    cache_.clear();
    sample_indices_.clear();
    unprime();
}

bool HaplotypeLikelihoodCache::is_primed() const noexcept
{
    return static_cast<bool>(primed_sample_);
}

void HaplotypeLikelihoodCache::prime(const SampleName& sample) const
{
    primed_sample_ = sample_indices_.at(sample);
}

void HaplotypeLikelihoodCache::unprime() const noexcept
{
    primed_sample_ = boost::none;
}

// private methods

void HaplotypeLikelihoodCache::set_read_iterators_and_sample_indices(const ReadMap& reads)
{
    read_iterators_.clear();
    sample_indices_.clear();
    const auto num_samples = reads.size();
    if (read_iterators_.capacity() < num_samples) {
        read_iterators_.reserve(num_samples);
    }
    if (sample_indices_.bucket_count() < num_samples) {
        sample_indices_.rehash(num_samples);
    }
    std::size_t i {0};
    for (const auto& p : reads) {
        read_iterators_.emplace_back(std::cbegin(p.second), std::cend(p.second));
        sample_indices_.emplace(p.first, i++);
    }
}

// non-member methods

HaplotypeLikelihoodCache merge_samples(const std::vector<SampleName>& samples,
                                       const SampleName& new_sample,
                                       const std::vector<Haplotype>& haplotypes,
                                       const HaplotypeLikelihoodCache& haplotype_likelihoods)
{
    HaplotypeLikelihoodCache result {static_cast<unsigned>(haplotypes.size()), {new_sample}};
    for (const auto& haplotype : haplotypes) {
        HaplotypeLikelihoodCache::LikelihoodVector likelihoods {};
        for (const auto& sample : samples) {
            const auto& m = haplotype_likelihoods(sample, haplotype);
            likelihoods.insert(std::end(likelihoods), std::cbegin(m), std::cend(m));
        }
        likelihoods.shrink_to_fit();
        result.insert(new_sample, haplotype, std::move(likelihoods));
    }
    return result;
}

namespace debug {

std::vector<std::reference_wrapper<const Haplotype>>
rank_haplotypes(const std::vector<Haplotype>& haplotypes, const SampleName& sample,
                const HaplotypeLikelihoodCache& haplotype_likelihoods)
{
    std::vector<std::pair<std::reference_wrapper<const Haplotype>, double>> ranks {};
    ranks.reserve(haplotypes.size());
    for (const auto& haplotype : haplotypes) {
        const auto& likelihoods = haplotype_likelihoods(sample, haplotype);
        ranks.emplace_back(haplotype, std::accumulate(std::cbegin(likelihoods), std::cend(likelihoods), 0.0));
    }
    std::sort(std::begin(ranks), std::end(ranks),
              [] (const auto& lhs, const auto& rhs) {
                  return lhs.second > rhs.second;
              });
    std::vector<std::reference_wrapper<const Haplotype>> result {};
    result.reserve(haplotypes.size());
    std::transform(std::cbegin(ranks), std::cend(ranks), std::back_inserter(result),
                   [] (const auto& p) { return p.first; });
    return result;
}

void print_read_haplotype_likelihoods(const std::vector<Haplotype>& haplotypes,
                                      const ReadMap& reads,
                                      const HaplotypeLikelihoodCache& haplotype_likelihoods,
                                      const std::size_t n)
{
    print_read_haplotype_likelihoods(std::cout, haplotypes, reads, haplotype_likelihoods, n);
}

} // namespace debug

} // namespace octopus
