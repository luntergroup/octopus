//
//  caching_fasta.cpp
//  Octopus
//
//  Created by Daniel Cooke on 22/07/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "caching_fasta.hpp"

#include <iterator>
#include <algorithm>
#include <stdexcept>
#include <utility>

#include "genomic_region.hpp"

#include <iostream> // TEST

// public methods

CachingFasta::CachingFasta(fs::path fasta_path)
:
fasta_ {std::move(fasta_path)},
contig_size_cache_ {},
sequence_cache_ {},
recently_used_regions_ {}
{
    setup_cache();
}

CachingFasta::CachingFasta(fs::path fasta_path, SizeType max_cache_size)
:
fasta_ {std::move(fasta_path)},
contig_size_cache_ {},
sequence_cache_ {},
recently_used_regions_ {},
max_cache_size_ {max_cache_size}
{
    setup_cache();
}

CachingFasta::CachingFasta(fs::path fasta_path, fs::path fasta_index_path)
:
fasta_ {std::move(fasta_path), std::move(fasta_index_path)},
contig_size_cache_ {},
sequence_cache_ {},
recently_used_regions_ {}
{
    setup_cache();
}

CachingFasta::CachingFasta(fs::path fasta_path, fs::path fasta_index_path, SizeType max_cache_size)
:
fasta_ {std::move(fasta_path), std::move(fasta_index_path)},
contig_size_cache_ {},
sequence_cache_ {},
recently_used_regions_ {},
max_cache_size_ {max_cache_size}
{
    setup_cache();
}

CachingFasta::CachingFasta(fs::path fasta_path, SizeType max_cache_size, double locality_bias, double forward_bias)
:
fasta_ {std::move(fasta_path)},
contig_size_cache_ {},
sequence_cache_ {},
recently_used_regions_ {},
max_cache_size_ {max_cache_size},
locality_bias_ {locality_bias},
forward_bias_ {forward_bias}
{
    if (locality_bias_ < 0 || locality_bias_ > 1) {
        throw std::domain_error {std::string("Invalid locality bias ") + std::to_string(locality_bias_) + std::string("; must be [0, 1]")};
    }
    
    if (forward_bias_ < 0 || forward_bias_ > 1) {
        throw std::domain_error {std::string("Invalid forward bias ") + std::to_string(forward_bias_) + std::string("; must be [0, 1]")};
    }
    
    setup_cache();
}

// virtual private methods

std::string CachingFasta::do_get_reference_name() const
{
    return fasta_.get_reference_name();
}

std::vector<std::string> CachingFasta::do_get_contig_names()
{
    return fasta_.get_contig_names();
}

CachingFasta::SizeType CachingFasta::do_get_contig_size(const std::string& contig_name)
{
    return contig_size_cache_.at(contig_name);
}

CachingFasta::SequenceType CachingFasta::do_fetch_sequence(const GenomicRegion& region)
{
    //std::cout << "requested " << region << std::endl;
    
    if (size(region) > max_cache_size_) {
        //std::cout << "region too big to cache.. fetching from file" << std::endl;
        return fasta_.fetch_sequence(region);
    } else if (is_contig_cached(region)) {
        auto it = get_cache_iterator(region);
        if (it != sequence_cache_.at(region.get_contig_name()).cend()) {
            //std::cout << "region is in cache contained by contig region " << it->first << std::endl;
            update_cache_position(region);
            return get_subsequence(region.get_contig_region(), it->first, it->second);
        }
    }
    
    auto fetch_region = get_region_to_fetch(region);
    
    //std::cout << region << " not in cache" << std::endl;
    //std::cout << "fetching " << fetch_region << " from file" << std::endl;
    
    auto fetched_sequence = fasta_.fetch_sequence(fetch_region);
    
    auto result = get_subsequence(region.get_contig_region(), fetch_region.get_contig_region(), fetched_sequence);
    
    add_sequence_to_cache(std::move(fetched_sequence), fetch_region);
    
    return result;
}

// non-virtual private methods

void CachingFasta::setup_cache()
{
    auto contig_names = fasta_.get_contig_names();
    
    contig_size_cache_.reserve(contig_names.size());
    
    for (auto&& contig_name : contig_names) {
        auto size = fasta_.get_contig_size(contig_name);
        genome_size_ += size;
        contig_size_cache_.emplace(std::move(contig_name), size);
    }
}

CachingFasta::SizeType CachingFasta::get_remaining_cache_size() const
{
    return max_cache_size_ - current_cache_size_;
}

bool CachingFasta::is_contig_cached(const GenomicRegion& region) const
{
    return sequence_cache_.count(region.get_contig_name()) != 0;
}

CachingFasta::CacheIterator CachingFasta::get_cache_iterator(const GenomicRegion& region) const
{
    // precondition that contig is in cache
    const auto& contig_cache = sequence_cache_.at(region.get_contig_name());
    
    if (contig_cache.empty()) return contig_cache.cend();
    
    const auto& contig_region = region.get_contig_region();
    
    auto it = contig_cache.lower_bound(contig_region);
    
    if (it == contig_cache.cbegin() || contig_region == it->first) {
        return contains(it->first, contig_region) ? it : contig_cache.cend();
    } else {
        return contains(std::prev(it)->first, contig_region) ? std::prev(it) : contig_cache.cend();
    }
}

GenomicRegion CachingFasta::get_region_to_fetch(const GenomicRegion& requested_region) const
{
    // We know the entire region is not in cache, but not parts of it may be
    if (sequence_cache_.count(requested_region.get_contig_name()) == 0) {
        return get_new_contig_chunk(requested_region);
    } else {
        return get_hit_contig_chunk(requested_region);
    }
}

CachingFasta::SizeType CachingFasta::get_lhs_extension_size(const GenomicRegion& requested_region) const
{
    return std::min(requested_region.get_begin(),
                    static_cast<SizeType>((max_cache_size_ - size(requested_region)) * locality_bias_ * (1.0 - forward_bias_)));
}

CachingFasta::SizeType CachingFasta::get_rhs_extension_size(const GenomicRegion& requested_region) const
{
    return std::min(contig_size_cache_.at(requested_region.get_contig_name()) - requested_region.get_end(),
                    static_cast<SizeType>((max_cache_size_ - size(requested_region)) * locality_bias_ * forward_bias_));
}

GenomicRegion CachingFasta::get_new_contig_chunk(const GenomicRegion& requested_region) const
{
    auto begin = requested_region.get_begin() - get_lhs_extension_size(requested_region);
    auto end   = requested_region.get_end() + get_rhs_extension_size(requested_region);
    
    return GenomicRegion {requested_region.get_contig_name(), begin, end};
}

GenomicRegion CachingFasta::get_hit_contig_chunk(const GenomicRegion& requested_region) const
{
    // TODO: optimise this to avoid requesting sequence we already have in cache
    
//    auto cached_overlapped = this->overlap_range(requested_region);
//    
//    auto num_overlapped = std::distance(cached_overlapped.begin(), cached_overlapped.end());
//    
//    if (num_overlapped == 0) {
//        std::cout << "no overlaps" << std::endl;
//    } else if (num_overlapped == 1) {
//        std::cout << "single overlap in cache" << std::endl;
//    } else {
//        std::cout << "multiple overlaps in cache" << std::endl;
//        // covers multiple fragements
//    }
    
    return get_new_contig_chunk(requested_region);
}

void CachingFasta::add_sequence_to_cache(SequenceType&& sequence, const GenomicRegion& region)
{
    recache_overlapped_regions(sequence, region);
    
    const auto& contig = region.get_contig_name();
    
    if (sequence_cache_.count(contig) == 0 && contig_size_cache_.at(contig) >= get_remaining_cache_size()) {
        auto target_cache_size = static_cast<SizeType>(current_cache_size_ * (1.0 - locality_bias_));
        
        while (current_cache_size_ > target_cache_size) {
            //std::cout << "optimizing cache... removing " << recently_used_regions_.back() << std::endl;
            remove_from_sequence_cache(recently_used_regions_.back());
            current_cache_size_ -= size(recently_used_regions_.back());
            recently_used_regions_.pop_back();
        }
    }
    
    sequence_cache_[region.get_contig_name()].emplace(region.get_contig_region(), std::move(sequence));
    recently_used_regions_.push_front(region);
    current_cache_size_ += size(region);
    
    while (current_cache_size_ > max_cache_size_) {
        //std::cout << "cache full... removing " << recently_used_regions_.back() << std::endl;
        remove_from_sequence_cache(recently_used_regions_.back());
        current_cache_size_ -= size(recently_used_regions_.back());
        recently_used_regions_.pop_back();
    }
    
    //std::cout << "added " << region << " to cache" << std::endl;
    //std::cout << "used cache size is " << current_cache_size_ << std::endl;
}

void CachingFasta::update_cache_position(const GenomicRegion& region)
{
    auto it = std::find_if(recently_used_regions_.cbegin(), recently_used_regions_.cend(),
                           [&region] (const auto& cached_region) {
                               return contains(cached_region, region);
                           });
    
    recently_used_regions_.push_front(*it);
    recently_used_regions_.erase(it);
}

CachingFasta::OverlapRange CachingFasta::overlap_range(const GenomicRegion& region) const
{
    // precondition that contig is in cache
    const auto& contig_cache = sequence_cache_.at(region.get_contig_name());
    
    const auto& contig_region = region.get_contig_region();
    
    auto it = contig_cache.lower_bound(contig_region);
    
    if (it != contig_cache.cbegin()) --it;
    
    return {it, std::find_if_not(it, contig_cache.cend(),
                                 [&contig_region] (const auto& cached_region) {
                                     return overlaps(contig_region, cached_region.first);
                                 })};
}

void CachingFasta::remove_from_sequence_cache(const GenomicRegion& region)
{
    //std::cout << "removing " << region << " from cache" << std::endl;
    sequence_cache_[region.get_contig_name()].erase(region.get_contig_region());
}

void CachingFasta::remove_from_usage_cache(const GenomicRegion& region)
{
    recently_used_regions_.remove(region);
}

void CachingFasta::replace_in_usage_cache(const GenomicRegion& from, const GenomicRegion& to)
{
    auto it = std::find(recently_used_regions_.begin(), recently_used_regions_.end(), from);
    recently_used_regions_.insert(recently_used_regions_.erase(it), to);
}

template <typename T>
T get_nonoverlapped(const T& lhs, const T& rhs)
{
    return (begins_before(lhs, rhs)) ? get_left_overhang(lhs, rhs) : get_right_overhang(lhs, rhs);
}

void CachingFasta::recache_overlapped_regions(const SequenceType& sequence, const GenomicRegion& region)
{
    if (sequence_cache_.count(region.get_contig_name()) == 1) {
        auto cached_overlapped = this->overlap_range(region);
        auto num_overlapped = std::distance(cached_overlapped.begin(), cached_overlapped.end());
        
        if (num_overlapped > 0) {
            std::vector<ContigSequenceCache::value_type> to_recache {};
            to_recache.reserve(num_overlapped);
            
            std::for_each(cached_overlapped.begin(), cached_overlapped.end(), [this, &region, &to_recache] (auto& p){
                const ContigRegion& overlapped_contig_region {p.first};
                const GenomicRegion overlapped_region {region.get_contig_name(), overlapped_contig_region};
                
                current_cache_size_ -= size(overlapped_region);
                
                if (contains(region.get_contig_region(), overlapped_contig_region)) {
                    remove_from_usage_cache(overlapped_region);
                    //std::cout << "removed " << overlapped_region << " from cache" << std::endl;
                } else {
                    auto nonoverlapped_region = get_nonoverlapped(overlapped_region, region);
                    const auto& nonoverlapped_contig_region = nonoverlapped_region.get_contig_region();
                    
                    to_recache.emplace_back(nonoverlapped_contig_region,
                                            get_subsequence(nonoverlapped_contig_region,
                                                            overlapped_contig_region, p.second));
                    
                    replace_in_usage_cache(overlapped_region, nonoverlapped_region);
                    current_cache_size_ += size(nonoverlapped_contig_region);
                    //std::cout << "replaced " << overlapped_region << " with " << nonoverlapped_region << " in cache" << std::endl;
                }
            });
            
            // update sequence_cache_ here as std::map iterators are invalidated on erase
            sequence_cache_.at(region.get_contig_name()).erase(cached_overlapped.begin(), cached_overlapped.end());
            sequence_cache_.at(region.get_contig_name()).insert(std::make_move_iterator(to_recache.begin()),
                                                                std::make_move_iterator(to_recache.end()));
        }
    }
}

CachingFasta::SequenceType CachingFasta::get_subsequence(const ContigRegion& requested_region,
                                                         const ContigRegion& sequence_region,
                                                         const SequenceType& sequence) const
{
    // precondition is contains(sequence_region, requested_region)
    return sequence.substr(requested_region.get_begin() - sequence_region.get_begin(), size(requested_region));
}
