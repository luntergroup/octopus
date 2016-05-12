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
#include <cassert>

#include "genomic_region.hpp"

#include <iostream> // TEST

// public methods

CachingFasta::CachingFasta(Path fasta_path)
:
fasta_ {std::move(fasta_path)},
contig_sizes_ {},
sequence_cache_ {},
recently_used_regions_ {}
{
    setup_cache();
}

CachingFasta::CachingFasta(Path fasta_path, SizeType max_cache_size)
:
fasta_ {std::move(fasta_path)},
contig_sizes_ {},
sequence_cache_ {},
recently_used_regions_ {},
max_cache_size_ {max_cache_size}
{
    setup_cache();
}

CachingFasta::CachingFasta(Path fasta_path, Path fasta_index_path)
:
fasta_ {std::move(fasta_path), std::move(fasta_index_path)},
contig_sizes_ {},
sequence_cache_ {},
recently_used_regions_ {}
{
    setup_cache();
}

CachingFasta::CachingFasta(Path fasta_path, Path fasta_index_path, SizeType max_cache_size)
:
fasta_ {std::move(fasta_path), std::move(fasta_index_path)},
contig_sizes_ {},
sequence_cache_ {},
recently_used_regions_ {},
max_cache_size_ {max_cache_size}
{
    setup_cache();
}

CachingFasta::CachingFasta(Path fasta_path, SizeType max_cache_size, double locality_bias, double forward_bias)
:
fasta_ {std::move(fasta_path)},
contig_sizes_ {},
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

bool CachingFasta::do_is_open() const noexcept
{
    return fasta_.is_open();
}

std::string CachingFasta::do_get_reference_name() const
{
    return fasta_.get_reference_name();
}

std::vector<CachingFasta::ContigNameType> CachingFasta::do_get_contig_names() const
{
    return fasta_.get_contig_names();
}

CachingFasta::SizeType CachingFasta::do_get_contig_size(const ContigNameType& contig) const
{
    return contig_sizes_.at(contig);
}

CachingFasta::SequenceType CachingFasta::do_fetch_sequence(const GenomicRegion& region) const
{
    if (size(region) > max_cache_size_) {
        return fasta_.fetch_sequence(region);
    } else if (is_contig_cached(region)) {
        const auto it = get_cache_iterator(region);
        if (it != std::cend(sequence_cache_.at(region.get_contig_name()))) {
            update_cache_position(region);
            return get_subsequence(region.get_contig_region(), it->first, it->second);
        }
    }
    
    auto fetch_region = get_region_to_fetch(region);
    
    auto fetched_sequence = fasta_.fetch_sequence(fetch_region);
    
    auto result = get_subsequence(region.get_contig_region(), fetch_region.get_contig_region(),
                                  fetched_sequence);
    
    add_sequence_to_cache(std::move(fetched_sequence), fetch_region);
    
    return result;
}

// non-virtual private methods

void CachingFasta::setup_cache()
{
    auto contig_names = fasta_.get_contig_names();
    
    contig_sizes_.reserve(contig_names.size());
    
    for (auto&& contig_name : contig_names) {
        auto size = fasta_.get_contig_size(contig_name);
        genome_size_ += size;
        contig_sizes_.emplace(std::move(contig_name), size);
    }
}

CachingFasta::SizeType CachingFasta::get_remaining_cache_size() const
{
    assert(max_cache_size_ >= current_cache_size_);
    return max_cache_size_ - current_cache_size_;
}

bool CachingFasta::is_contig_cached(const GenomicRegion& region) const
{
    return sequence_cache_.count(region.get_contig_name()) != 0;
}

CachingFasta::CacheIterator CachingFasta::get_cache_iterator(const GenomicRegion& region) const
{
    assert(sequence_cache_.count(region.get_contig_name()) == 1);
    
    const auto& contig_cache = sequence_cache_.at(region.get_contig_name());
    
    if (contig_cache.empty()) return std::cend(contig_cache);
    
    const auto& contig_region = region.get_contig_region();
    
    const auto it = contig_cache.lower_bound(contig_region);
    
    if (it == std::cend(contig_cache)) return it;
    
    if (it == std::cbegin(contig_cache) || contig_region == it->first) {
        return contains(it->first, contig_region) ? it : std::cend(contig_cache);
    } else {
        return contains(std::prev(it)->first, contig_region) ? std::prev(it) : std::cend(contig_cache);
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
    return std::min(contig_sizes_.at(requested_region.get_contig_name()) - requested_region.get_end(),
                    static_cast<SizeType>((max_cache_size_ - size(requested_region)) * locality_bias_ * forward_bias_));
}

GenomicRegion CachingFasta::get_new_contig_chunk(const GenomicRegion& requested_region) const
{
    const auto begin = requested_region.get_begin() - get_lhs_extension_size(requested_region);
    const auto end   = requested_region.get_end() + get_rhs_extension_size(requested_region);
    
    return GenomicRegion {requested_region.get_contig_name(), begin, end};
}

GenomicRegion CachingFasta::get_hit_contig_chunk(const GenomicRegion& requested_region) const
{
    // TODO: optimise this to avoid requesting sequence we already have in cache
    return get_new_contig_chunk(requested_region);
}

void CachingFasta::add_sequence_to_cache(SequenceType&& sequence, const GenomicRegion& region) const
{
    recache_overlapped_regions(sequence, region);
    
    const auto& contig = region.get_contig_name();
    
    if (sequence_cache_.count(contig) == 0 && contig_sizes_.at(contig) >= get_remaining_cache_size()) {
        const auto target_cache_size = static_cast<SizeType>(current_cache_size_ * (1.0 - locality_bias_));
        
        while (current_cache_size_ > target_cache_size) {
            remove_from_sequence_cache(recently_used_regions_.back());
            current_cache_size_ -= size(recently_used_regions_.back());
            recently_used_regions_.pop_back();
        }
    }
    
    sequence_cache_[region.get_contig_name()].emplace(region.get_contig_region(), std::move(sequence));
    recently_used_regions_.push_front(region);
    current_cache_size_ += size(region);
    
    while (current_cache_size_ > max_cache_size_) {
        remove_from_sequence_cache(recently_used_regions_.back());
        current_cache_size_ -= size(recently_used_regions_.back());
        recently_used_regions_.pop_back();
    }
}

void CachingFasta::update_cache_position(const GenomicRegion& region) const
{
    const auto it = std::find_if(std::cbegin(recently_used_regions_), std::cend(recently_used_regions_),
                                 [&region] (const auto& cached_region) {
                                     return contains(cached_region, region);
                                 });
    
    recently_used_regions_.push_front(*it);
    recently_used_regions_.erase(it);
}

CachingFasta::OverlapRange CachingFasta::overlap_range(const GenomicRegion& region) const
{
    assert(sequence_cache_.count(region.get_contig_name()) == 1);
    
    const auto& contig_cache = sequence_cache_.at(region.get_contig_name());
    
    const auto& contig_region = region.get_contig_region();
    
    auto it = contig_cache.lower_bound(contig_region);
    
    if (it != std::cbegin(contig_cache)) --it;
    
    return {it, std::find_if_not(it, std::cend(contig_cache),
                                 [&contig_region] (const auto& cached_region) {
                                     return overlaps(contig_region, cached_region.first);
                                 })};
}

void CachingFasta::remove_from_sequence_cache(const GenomicRegion& region) const
{
    sequence_cache_[region.get_contig_name()].erase(region.get_contig_region());
}

void CachingFasta::remove_from_usage_cache(const GenomicRegion& region) const
{
    recently_used_regions_.remove(region);
}

void CachingFasta::replace_in_usage_cache(const GenomicRegion& from, const GenomicRegion& to) const
{
    const auto it = std::find(std::begin(recently_used_regions_), std::end(recently_used_regions_), from);
    recently_used_regions_.insert(recently_used_regions_.erase(it), to);
}

template <typename T>
T get_nonoverlapped(const T& lhs, const T& rhs)
{
    return (begins_before(lhs, rhs)) ? left_overhang_region(lhs, rhs) : right_overhang_region(lhs, rhs);
}

void CachingFasta::recache_overlapped_regions(const SequenceType& sequence, const GenomicRegion& region) const
{
    if (sequence_cache_.count(region.get_contig_name()) == 1) {
        auto cached_overlapped = this->overlap_range(region);
        
        const auto num_overlapped = std::distance(std::cbegin(cached_overlapped), std::cend(cached_overlapped));
        
        if (num_overlapped > 0) {
            std::vector<ContigSequenceCache::value_type> to_recache {};
            to_recache.reserve(num_overlapped);
            
            std::for_each(std::begin(cached_overlapped), std::end(cached_overlapped),
                          [this, &region, &to_recache] (auto& p){
                const ContigRegion& overlapped_contig_region {p.first};
                const GenomicRegion overlapped_region {region.get_contig_name(), overlapped_contig_region};
                
                current_cache_size_ -= size(overlapped_region);
                
                if (contains(region.get_contig_region(), overlapped_contig_region)) {
                    remove_from_usage_cache(overlapped_region);
                } else {
                    auto nonoverlapped_region = get_nonoverlapped(overlapped_region, region);
                    const auto& nonoverlapped_contig_region = nonoverlapped_region.get_contig_region();
                    
                    to_recache.emplace_back(nonoverlapped_contig_region,
                                            get_subsequence(nonoverlapped_contig_region,
                                                            overlapped_contig_region, p.second));
                    
                    replace_in_usage_cache(overlapped_region, nonoverlapped_region);
                    current_cache_size_ += size(nonoverlapped_contig_region);
                }
            });
            
            // update sequence_cache_ here as std::map iterators are invalidated on erase
            sequence_cache_.at(region.get_contig_name()).erase(std::begin(cached_overlapped),
                                                               std::end(cached_overlapped));
            sequence_cache_.at(region.get_contig_name()).insert(std::make_move_iterator(std::begin(to_recache)),
                                                                std::make_move_iterator(std::end(to_recache)));
        }
    }
}

CachingFasta::SequenceType CachingFasta::get_subsequence(const ContigRegion& requested_region,
                                                         const ContigRegion& sequence_region,
                                                         const SequenceType& sequence) const
{
    assert(contains(sequence_region, requested_region));
    return sequence.substr(begin_distance(requested_region, sequence_region), size(requested_region));
}
