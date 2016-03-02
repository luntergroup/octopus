//
//  haplotype_likelihood_cache.cpp
//  Octopus
//
//  Created by Daniel Cooke on 25/11/2015.
//  Copyright © 2015 Oxford University. All rights reserved.
//

#include "haplotype_likelihood_cache.hpp"

#include <vector>
#include <algorithm>
#include <iostream>
#include <utility>

#include <boost/iterator/transform_iterator.hpp>

namespace Octopus
{
// public methods

HaplotypeLikelihoodCache::HaplotypeLikelihoodCache(const unsigned max_haplotypes,
                                                   const std::vector<SampleIdType>& samples)
:
cache_ {max_haplotypes},
sample_indices_ {samples.size()}
{}

void HaplotypeLikelihoodCache::populate(const ReadMap& reads,
                                        const std::vector<Haplotype>& haplotypes,
                                        HaplotypeLikelihoodModel::InactiveRegionState flank_state)
{
    if (!cache_.empty()) {
        cache_.clear();
    }
    
    if (cache_.bucket_count() < haplotypes.size()) {
        cache_.rehash(haplotypes.size());
    }
    
    sample_indices_.clear();
    
    set_read_iterators_and_sample_indices(reads);
    
    const auto num_samples = reads.size();
    
    std::vector<std::vector<KmerPerfectHashes>> read_hashes {};
    read_hashes.reserve(num_samples);
    
    for (const auto& t : read_iterators_) {
        std::vector<KmerPerfectHashes> sample_read_hashes {};
        sample_read_hashes.reserve(t.num_reads);
        
        std::transform(t.first, t.last, std::back_inserter(sample_read_hashes),
                       [] (const AlignedRead& read) {
                           return compute_kmer_hashes<KMER_SIZE>(read.get_sequence());
                       });
        
        read_hashes.emplace_back(std::move(sample_read_hashes));
    }
    
    auto haplotype_hashes = init_kmer_hash_table<KMER_SIZE>();
    
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
        populate_kmer_hash_table<KMER_SIZE>(haplotype.get_sequence(), haplotype_hashes);
        
        auto haplotype_mapping_counts = init_mapping_counts(haplotype_hashes);
        
        auto it = std::begin(cache_.emplace(std::piecewise_construct,
                                            std::forward_as_tuple(haplotype),
                                            std::forward_as_tuple(num_samples)).first->second);
        
        const auto read_hash_itr = std::cbegin(read_hashes);
        
        HaplotypeLikelihoodModel likelihood_model {haplotype, flank_state};
        
        for (const auto& t : read_iterators_) {
            *it = std::vector<double>(t.num_reads);
            
            std::transform(t.first, t.last, std::cbegin(*read_hash_itr), std::begin(*it),
                           [&] (const AlignedRead& read, const auto& read_hashes) {
                               const auto last_mapping_position = map_query_to_target(read_hashes, haplotype_hashes,
                                                                                      haplotype_mapping_counts,
                                                                                      first_mapping_position);
                               
                               reset_mapping_counts(haplotype_mapping_counts);
                               
                               return likelihood_model.log_probability(read, first_mapping_position,
                                                                       last_mapping_position);
                           });
            
            ++it;
        }
        
        clear_kmer_hash_table(haplotype_hashes);
    }
    
    read_iterators_.clear();
}

const HaplotypeLikelihoodCache::ReadProbabilities&
HaplotypeLikelihoodCache::log_likelihoods(const SampleIdType& sample,
                                          const Haplotype& haplotype) const
{
    return cache_.at(haplotype)[sample_indices_.at(sample)];
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
    
//    void HaplotypeLikelihoodCache::set_read_hashes()
//    {
//        if (read_hashes_.size() < sample_indices_.size()) {
//            read_hashes_.resize(sample_indices_.size());
//        }
//        
//        std::transform(std::cbegin(read_iterators_), std::cend(read_iterators_),
//                       std::begin(read_hashes_),
//                       [] (const auto& t) {
//                           
//                       });
//        
//        for (const auto& t : read_iterators_) {
//            std::vector<KmerPerfectHashes> sample_read_hashes {};
//            sample_read_hashes.reserve(t.num_reads);
//            
//            std::transform(t.first, t.last, std::back_inserter(sample_read_hashes),
//                           [] (const AlignedRead& read) {
//                               return compute_kmer_hashes<KMER_SIZE>(read.get_sequence());
//                           });
//            
//            read_hashes_.emplace_back(std::move(sample_read_hashes));
//        }
//    }

namespace debug
{
//    void print_read_haplotype_liklihoods(const std::vector<Haplotype>& haplotypes,
//                                         const ReadMap& reads,
//                                         HaplotypeLikelihoodCache& haplotype_likelihoods,
//                                         size_t n)
//    {
//        auto m = std::min(n, haplotypes.size());
//        
//        std::cout << "debug: top " << m << " haplotype likelihoods for each read in each sample" << std::endl;
//        
//        using HaplotypeReference = std::reference_wrapper<const Haplotype>;
//        
//        for (const auto& sample_reads : reads) {
//            std::cout << "Sample: " << sample_reads.first << ":" << std::endl;
//            
//            for (const auto& read : sample_reads.second) {
//                std::cout << "\tRead: " << read.get_region() << " " << read.get_cigar_string() << ":" << std::endl;
//                
//                std::vector<std::pair<HaplotypeReference, double>> top {};
//                top.reserve(haplotypes.size());
//                
//                for (const auto& haplotype : haplotypes) {
//                    top.emplace_back(haplotype, haplotype_likelihoods.log_probability(read, haplotype));
//                }
//                
//                std::sort(std::begin(top), std::end(top),
//                          [] (const auto& lhs, const auto& rhs) {
//                              return lhs.second > rhs.second;
//                          });
//                
//                for (unsigned i {}; i < m; ++i) {
//                    std::cout << "\t\t ";
//                    print_variant_alleles(top[i].first);
//                    std::cout << " " << std::setprecision(10) << top[i].second << std::endl;
//                }
//            }
//        }
//    }
} // namespace debug
} // namespace Octopus
