//
//  haplotype_likelihood_cache.cpp
//  Octopus
//
//  Created by Daniel Cooke on 25/11/2015.
//  Copyright © 2015 Oxford University. All rights reserved.
//

#include "haplotype_likelihood_cache.hpp"

#include <utility>
#include <cassert>

#include <iostream> // DEBUG
#include <iomanip>  // DEBUG

namespace Octopus
{
// public methods

HaplotypeLikelihoodCache::HaplotypeLikelihoodCache(const unsigned max_haplotypes,
                                                   const std::vector<SampleIdType>& samples)
:
cache_ {max_haplotypes},
sample_indices_ {samples.size()}
{}

HaplotypeLikelihoodCache::HaplotypeLikelihoodCache(HaplotypeLikelihoodModel likelihood_model,
                                                   unsigned max_haplotypes,
                                                   const std::vector<SampleIdType>& samples)
:
likelihood_model_ {std::move(likelihood_model)},
cache_ {max_haplotypes},
sample_indices_ {samples.size()}
{}

HaplotypeLikelihoodCache::ReadPacket::ReadPacket(Iterator first, Iterator last)
:
first {first},
last {last},
num_reads {static_cast<std::size_t>(std::distance(first, last))}
{}

void HaplotypeLikelihoodCache::populate(const ReadMap& reads,
                                        const std::vector<Haplotype>& haplotypes,
                                        boost::optional<FlankState> flank_state)
{
    // This code is not very pretty because it is a real bottleneck for the entire application.
    // We want to try a minimise memory allocations for the mapping.
    
    if (!cache_.empty()) {
        cache_.clear();
    }
    
    if (cache_.bucket_count() < haplotypes.size()) {
        cache_.rehash(haplotypes.size());
    }
    
    sample_indices_.clear();
    
    set_read_iterators_and_sample_indices(reads);
    
    assert(reads.size() == read_iterators_.size());
    
    const auto num_samples = reads.size();
    
    std::vector<std::vector<KmerPerfectHashes>> read_hashes {};
    read_hashes.reserve(num_samples);
    
    for (const auto& t : read_iterators_) {
        std::vector<KmerPerfectHashes> sample_read_hashes {};
        sample_read_hashes.reserve(t.num_reads);
        
        std::transform(t.first, t.last, std::back_inserter(sample_read_hashes),
                       [] (const AlignedRead& read) {
                           return compute_kmer_hashes<MAPPER_KMER_SIZE>(read.get_sequence());
                       });
        
        read_hashes.emplace_back(std::move(sample_read_hashes));
    }
    
    auto haplotype_hashes = init_kmer_hash_table<MAPPER_KMER_SIZE>();
    
    const auto max_mapping_positions = sequence_size(*std::max_element(std::cbegin(haplotypes),
                                                                       std::cend(haplotypes),
                                                                       [] (const auto& lhs, const auto& rhs) {
                                                                           return sequence_size(lhs) < sequence_size(rhs);
                                                                       }));
    
    if (max_mapping_positions > mapping_positions_.size()) {
        mapping_positions_.resize(max_mapping_positions);
    }
    
    const auto first_mapping_position = std::begin(mapping_positions_);
    
    for (const auto& haplotype : haplotypes) {
        populate_kmer_hash_table<MAPPER_KMER_SIZE>(haplotype.get_sequence(), haplotype_hashes);
        
        auto haplotype_mapping_counts = init_mapping_counts(haplotype_hashes);
        
        auto it = std::begin(cache_.emplace(std::piecewise_construct,
                                            std::forward_as_tuple(haplotype),
                                            std::forward_as_tuple(num_samples)).first->second);
        
        auto read_hash_itr = std::cbegin(read_hashes);
        
        likelihood_model_.buffer(haplotype, flank_state);
        
        for (const auto& t : read_iterators_) {
            *it = std::vector<double>(t.num_reads);
            
            std::transform(t.first, t.last, std::cbegin(*read_hash_itr), std::begin(*it),
                           [&] (const AlignedRead& read, const auto& read_hashes) {
                               const auto last_mapping_position = map_query_to_target(read_hashes, haplotype_hashes,
                                                                                      haplotype_mapping_counts,
                                                                                      first_mapping_position);
                               
                               reset_mapping_counts(haplotype_mapping_counts);
                               
                               return likelihood_model_.log_probability(read, first_mapping_position,
                                                                        last_mapping_position);
                           });
            
            ++read_hash_itr;
            ++it;
        }
        
        clear_kmer_hash_table(haplotype_hashes);
    }
    
    likelihood_model_.clear();
    read_iterators_.clear();
}

const HaplotypeLikelihoodCache::Likelihoods&
HaplotypeLikelihoodCache::log_likelihoods(const SampleIdType& sample,
                                          const Haplotype& haplotype) const
{
    return cache_.at(haplotype)[sample_indices_.at(sample)];
}

HaplotypeLikelihoodCache::SampleLikelihoodMap
HaplotypeLikelihoodCache::extract_sample(const SampleIdType& sample) const
{
    const auto sample_index = sample_indices_.at(sample);
    
    SampleLikelihoodMap result {cache_.size()};
    
    for (const auto& p : cache_) {
        result.emplace(std::cref(p.first), std::cref(p.second[sample_index]));
    }
    
    return result;
}

void HaplotypeLikelihoodCache::reserve(std::size_t num_samples, std::size_t num_haplotypes)
{
    sample_indices_.reserve(num_samples);
    cache_.reserve(num_haplotypes);
}

void HaplotypeLikelihoodCache::clear()
{
    cache_.clear();
}

// private methods

void HaplotypeLikelihoodCache::set_read_iterators_and_sample_indices(const ReadMap& reads)
{
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

HaplotypeLikelihoodCache merge_samples(const std::vector<SampleIdType>& samples,
                                       const SampleIdType& new_sample,
                                       const std::vector<Haplotype>& haplotypes,
                                       const HaplotypeLikelihoodCache& haplotype_likelihoods)
{
    HaplotypeLikelihoodCache result {};
    result.reserve(1, haplotypes.size());
    
    for (const auto& haplotype : haplotypes) {
        HaplotypeLikelihoodCache::Likelihoods likelihoods {};
        for (const auto& sample : samples) {
            const auto& m = haplotype_likelihoods.log_likelihoods(sample, haplotype);
            likelihoods.insert(std::end(likelihoods), std::cbegin(m), std::cend(m));
        }
        likelihoods.shrink_to_fit();
        result.insert(new_sample, haplotype, std::move(likelihoods));
    }
    
    return result;
}

namespace debug
{
    std::vector<std::reference_wrapper<const Haplotype>>
    rank_haplotypes(const std::vector<Haplotype>& haplotypes, const SampleIdType& sample,
                    const HaplotypeLikelihoodCache& haplotype_likelihoods)
    {
        std::vector<std::pair<std::reference_wrapper<const Haplotype>, double>> ranks {};
        ranks.reserve(haplotypes.size());
        
        for (const auto& haplotype : haplotypes) {
            const auto& likelihoods = haplotype_likelihoods.log_likelihoods(sample, haplotype);
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
    
    void print_read_haplotype_liklihoods(const std::vector<Haplotype>& haplotypes,
                                         const ReadMap& reads,
                                         const HaplotypeLikelihoodCache& haplotype_likelihoods,
                                         const std::size_t n)
    {
        print_read_haplotype_liklihoods(std::cout, haplotypes, reads, haplotype_likelihoods);
    }
} // namespace debug
} // namespace Octopus
