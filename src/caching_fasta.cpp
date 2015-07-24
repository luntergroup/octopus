//
//  caching_fasta.cpp
//  Octopus
//
//  Created by Daniel Cooke on 22/07/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "caching_fasta.h"

#include <iterator>
#include <algorithm>

#include "genomic_region.h"

#include <iostream> // TEST

CachingFasta::CachingFasta(std::string fasta_path)
:
fasta_ {std::move(fasta_path)},
contig_size_cache_ {},
sequence_cache_ {},
recently_used_regions_ {},
used_cache_size_ {},
max_cache_size_ {1000000}
{
    setup_cache();
}

CachingFasta::CachingFasta(std::string fasta_path, SizeType max_cache_size)
:
fasta_ {std::move(fasta_path)},
contig_size_cache_ {},
sequence_cache_ {},
recently_used_regions_ {},
used_cache_size_ {},
max_cache_size_ {max_cache_size}
{
    setup_cache();
}

CachingFasta::CachingFasta(std::string fasta_path, std::string fasta_index_path)
:
fasta_ {std::move(fasta_path), std::move(fasta_index_path)},
contig_size_cache_ {},
sequence_cache_ {},
recently_used_regions_ {},
used_cache_size_ {},
max_cache_size_ {1000000}
{
    setup_cache();
}

CachingFasta::CachingFasta(std::string fasta_path, std::string fasta_index_path, SizeType max_cache_size)
:
fasta_ {std::move(fasta_path), std::move(fasta_index_path)},
contig_size_cache_ {},
sequence_cache_ {},
recently_used_regions_ {},
used_cache_size_ {},
max_cache_size_ {max_cache_size}
{
    setup_cache();
}

std::string CachingFasta::get_reference_name()
{
    return fasta_.get_reference_name();
}

std::vector<std::string> CachingFasta::get_contig_names()
{
    return fasta_.get_contig_names();
}

CachingFasta::SizeType CachingFasta::get_contig_size(const std::string& contig_name)
{
    return contig_size_cache_.at(contig_name);
}

CachingFasta::SequenceType CachingFasta::get_sequence(const GenomicRegion& region)
{
    std::cout << "requested " << region << std::endl;
    
    if (size(region) > max_cache_size_) {
        std::cout << "too big to cache" << std::endl;
        return fasta_.get_sequence(region);
    } else if (is_region_cached(region)) {
        std::cout << "is in cache" << std::endl;
        update_cache_position(region);
    } else {
        std::cout << "not in cache" << std::endl;
        auto fetch_region = region_to_fetch(region);
        //std::cout << "fetching " << fetch_region << " from file" << std::endl;
        add_sequence_to_cache(fasta_.get_sequence(fetch_region), fetch_region);
    }
    
    CacheIterator it {get_cache_iterator(region)};
    
    return get_subsequence(region, it->first, it->second);
}

// private methods

void CachingFasta::setup_cache()
{
    for (const auto& contig_name : fasta_.get_contig_names()) {
        contig_size_cache_.emplace(contig_name, fasta_.get_contig_size(contig_name));
    }
}

// TODO: we should also check if we can chain together adjacent cached regions to make
// the requested region. This could save a lot of recaching.
bool CachingFasta::is_region_cached(const GenomicRegion& region) const
{
    if (sequence_cache_.count(region.get_contig_name()) == 0) return false;
    
    const auto& contig_cache = sequence_cache_.at(region.get_contig_name());
    
    auto it = contig_cache.lower_bound(region);
    
    if (it == contig_cache.cbegin() || region == it->first) {
        return contains(it->first, region);
    } else {
        return contains(std::prev(it)->first, region);
    }
}

GenomicRegion CachingFasta::region_to_fetch(const GenomicRegion& requested_region) const
{
    auto remaining_cache_size = max_cache_size_ - used_cache_size_;
    auto requested_size       = size(requested_region);
    
    if (remaining_cache_size <= requested_size) return requested_region;
    
    //auto num_contigs          = contig_size_cache_.size();
    auto contig_size = contig_size_cache_.at(requested_region.get_contig_name());
    auto lhs_size    = requested_region.get_begin();
    auto rhs_size    = contig_size - requested_region.get_end();
    
    if (sequence_cache_.count(requested_region.get_contig_name()) == 0) {
        auto balance = std::min(remaining_cache_size - requested_size, static_cast<SizeType>(10000));
        
        while (balance > 0 && (lhs_size > 0 || rhs_size > 0)) {
            if (lhs_size > 0) { --lhs_size; --balance; }
            if (rhs_size > 0) { --rhs_size; --balance; }
        }
        
        return GenomicRegion {requested_region.get_contig_name(), lhs_size, contig_size - rhs_size};
    } else {
        // need to look what's already in cache
    }
    
    return requested_region;
}

void CachingFasta::add_sequence_to_cache(const SequenceType& sequence, const GenomicRegion& region)
{
    recache_overlapped_regions(sequence, region);
    
    sequence_cache_[region.get_contig_name()].emplace(region, sequence);
    recently_used_regions_.push_front(region);
    used_cache_size_ += size(region);
    
    while (used_cache_size_ > max_cache_size_) {
        std::cout << "cache full... removing " << recently_used_regions_.back() << std::endl;
        remove_from_cache(recently_used_regions_.back());
        used_cache_size_ -= size(recently_used_regions_.back());
        recently_used_regions_.pop_back();
    }
    
    std::cout << "added " << region << " to cache" << std::endl;
    std::cout << "used cache size is " << used_cache_size_ << std::endl;
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

CachingFasta::CacheIterator CachingFasta::get_cache_iterator(const GenomicRegion& region) const
{
    const auto& contig_cache = sequence_cache_.at(region.get_contig_name());
    
    auto it = contig_cache.lower_bound(region);
    
    return (it == contig_cache.cbegin() || region == it->first) ? it : std::prev(it);
}

std::pair<CachingFasta::CacheIterator, CachingFasta::CacheIterator> CachingFasta::overlap_range(const GenomicRegion& region) const
{
    const auto& contig_cache = sequence_cache_.at(region.get_contig_name());
    
    auto it = contig_cache.lower_bound(region);
    
    if (it != contig_cache.cbegin()) --it;
    
    return {it, std::find_if_not(it, contig_cache.cend(),
                                 [&region] (const auto& cached_region) {
                                     return overlaps(region, cached_region.first);
                                 })};
}

void CachingFasta::remove_from_cache(const GenomicRegion& region)
{
    //std::cout << "removing " << region << " from cache" << std::endl;
    sequence_cache_[region.get_contig_name()].erase(region);
}

void CachingFasta::recache_overlapped_regions(const SequenceType& sequence, const GenomicRegion& region)
{
    if (sequence_cache_.count(region.get_contig_name()) == 1) {
        auto cached_overlapped = this->overlap_range(region);
        auto num_overlapped = std::distance(cached_overlapped.first, cached_overlapped.second);
        
        if (num_overlapped > 0) {
            std::vector<std::pair<GenomicRegion, SequenceType>> to_recache {};
            to_recache.reserve(num_overlapped);
            
            std::for_each(cached_overlapped.first, cached_overlapped.second, [this, &region, &to_recache] (auto& p){
                const GenomicRegion& overlapped {p.first};
                used_cache_size_ -= size(overlapped);
                
                if (contains(region, overlapped)) {
                    recently_used_regions_.remove(overlapped);
                    //std::cout << "removed " << overlapped << " from cache" << std::endl;
                } else {
                    auto nonoverlapped_region = (overlapped < region) ? get_left_overhang(overlapped, region) : get_right_overhang(overlapped, region);
                    to_recache.emplace_back(nonoverlapped_region, get_subsequence(nonoverlapped_region, overlapped, p.second));
                    used_cache_size_ += size(nonoverlapped_region);
                    recently_used_regions_.insert(recently_used_regions_.erase(std::find(recently_used_regions_.begin(),
                                                                                         recently_used_regions_.end(), overlapped)),
                                                  nonoverlapped_region);
                    //std::cout << "replaced " << overlapped << " with " << nonoverlapped_region << " in cache" << std::endl;
                }
            });
            
            // update sequence_cache_ here as std::map iterators are invalidated on erase
            sequence_cache_.at(region.get_contig_name()).erase(cached_overlapped.first, cached_overlapped.second);
            sequence_cache_.at(region.get_contig_name()).insert(std::make_move_iterator(to_recache.begin()),
                                                                std::make_move_iterator(to_recache.end()));
        }
    }
}

CachingFasta::SequenceType CachingFasta::get_subsequence(const GenomicRegion& requested_region,
                                                         const GenomicRegion& sequence_region,
                                                         const SequenceType& sequence) const
{
    // we know contains(sequence_region, requested_region) is true
    auto offset = requested_region.get_begin() - sequence_region.get_begin();
    return sequence.substr(offset, size(requested_region));
}
