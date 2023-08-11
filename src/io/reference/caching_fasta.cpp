// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "caching_fasta.hpp"

#include <iterator>
#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <utility>
#include <cassert>

#include "basics/genomic_region.hpp"

namespace octopus { namespace io {

// public methods

CachingFasta::CachingFasta(std::unique_ptr<ReferenceReader> fasta)
: CachingFasta {std::move(fasta), 0}
{
    max_cache_size_ = std::accumulate(std::cbegin(contig_sizes_), std::cend(contig_sizes_), 0,
                                      [] (const auto curr, const auto& p) { return curr + p.second; });
}

CachingFasta::CachingFasta(std::unique_ptr<ReferenceReader> fasta,
                           GenomicSize max_cache_size)
: CachingFasta {std::move(fasta), max_cache_size, 0.99, 0.5}
{}

CachingFasta::CachingFasta(std::unique_ptr<ReferenceReader> fasta,
                           const GenomicSize max_cache_size,
                           const double locality_bias,
                           const double forward_bias)
: fasta_ {std::move(fasta)}
, contig_sizes_ {}
, sequence_cache_ {}
, recently_used_regions_ {}
, genome_size_ {0}
, max_cache_size_ {max_cache_size}
, current_cache_size_ {0}
, locality_bias_ {locality_bias}
, forward_bias_ {forward_bias}
{
    if (locality_bias_ < 0 || locality_bias_ > 1) {
        throw std::domain_error {std::string("Invalid locality bias ")
                + std::to_string(locality_bias_) + std::string("; must be [0, 1]")};
    }
    
    if (forward_bias_ < 0 || forward_bias_ > 1) {
        throw std::domain_error {std::string("Invalid forward bias ")
                + std::to_string(forward_bias_) + std::string("; must be [0, 1]")};
    }
    
    setup_cache();
}

CachingFasta::CachingFasta(const CachingFasta& other)
{
    std::lock_guard<std::mutex> lock {other.mutex_};
    fasta_              = other.fasta_->clone();
    contig_sizes_       = other.contig_sizes_;
    genome_size_        = other.genome_size_;
    max_cache_size_     = other.max_cache_size_;
    current_cache_size_ = other.current_cache_size_;
    locality_bias_      = other.locality_bias_;
    forward_bias_       = other.forward_bias_;
}

CachingFasta& CachingFasta::operator=(CachingFasta other)
{
    using std::swap;
    swap(fasta_             , other.fasta_);
    swap(contig_sizes_      , other.contig_sizes_);
    swap(genome_size_       , other.genome_size_);
    swap(max_cache_size_    , other.max_cache_size_);
    swap(current_cache_size_, other.current_cache_size_);
    swap(locality_bias_     , other.locality_bias_);
    swap(forward_bias_      , other.forward_bias_);
    return *this;
}

CachingFasta::CachingFasta(CachingFasta&& other)
{
    std::lock_guard<std::mutex> lock {other.mutex_};
    fasta_              = std::move(other.fasta_);
    contig_sizes_       = std::move(other.contig_sizes_);
    genome_size_        = std::move(other.genome_size_);
    max_cache_size_     = std::move(other.max_cache_size_);
    current_cache_size_ = std::move(other.current_cache_size_);
    locality_bias_      = std::move(other.locality_bias_);
    forward_bias_       = std::move(other.forward_bias_);
}

CachingFasta& CachingFasta::operator=(CachingFasta&& other)
{
    if (this != &other) {
        std::unique_lock<std::mutex> lock_lhs {mutex_, std::defer_lock}, lock_rhs {other.mutex_, std::defer_lock};
        std::lock(lock_lhs, lock_rhs);
        fasta_              = std::move(other.fasta_);
        contig_sizes_       = std::move(other.contig_sizes_);
        genome_size_        = std::move(other.genome_size_);
        max_cache_size_     = std::move(other.max_cache_size_);
        current_cache_size_ = std::move(other.current_cache_size_);
        locality_bias_      = std::move(other.locality_bias_);
        forward_bias_       = std::move(other.forward_bias_);
    }
    return *this;
}

// virtual private methods

std::unique_ptr<ReferenceReader> CachingFasta::do_clone() const
{
    return std::make_unique<CachingFasta>(*this);
}

bool CachingFasta::do_is_open() const noexcept
{
    return fasta_->is_open();
}

std::string CachingFasta::do_fetch_reference_name() const
{
    return fasta_->fetch_reference_name();
}

std::vector<CachingFasta::ContigName> CachingFasta::do_fetch_contig_names() const
{
    return fasta_->fetch_contig_names();
}

CachingFasta::GenomicSize CachingFasta::do_fetch_contig_size(const ContigName& contig) const
{
    return contig_sizes_.at(contig);
}

CachingFasta::GeneticSequence CachingFasta::do_fetch_sequence(const GenomicRegion& region) const
{
    if (is_empty(region)) {
        return "";
    }
    if (size(region) > max_cache_size_) {
        return fasta_->fetch_sequence(region);
    }
    std::unique_lock<std::mutex> lock {mutex_};
    const auto cache_itr = find_cached(region);
    if (cache_itr) {
        register_cache_hit(region);
        return get_subsequence(region.contig_region(), (*cache_itr)->first, (*cache_itr)->second);
    }
    auto fetch_region = get_region_to_fetch(region);
    assert(contains(fetch_region, region));
    lock.unlock();
    auto fetched_sequence = fasta_->fetch_sequence(fetch_region);
    auto result = get_subsequence(region.contig_region(), fetch_region.contig_region(), fetched_sequence);
    lock.lock();
    add_sequence_to_cache(std::move(fetched_sequence), std::move(fetch_region));
    return result;
}

// non-virtual private methods

void CachingFasta::setup_cache()
{
    auto contig_names = fasta_->fetch_contig_names();
    contig_sizes_.reserve(contig_names.size());
    for (auto&& contig_name : contig_names) {
        const auto size = fasta_->fetch_contig_size(contig_name);
        genome_size_ += size;
        contig_sizes_.emplace(std::move(contig_name), size);
    }
}

CachingFasta::GenomicSize CachingFasta::get_remaining_cache_size() const
{
    assert(max_cache_size_ >= current_cache_size_);
    return max_cache_size_ - current_cache_size_;
}

bool CachingFasta::is_contig_cached(const GenomicRegion& region) const
{
    return sequence_cache_.count(region.contig_name()) != 0;
}

boost::optional<CachingFasta::CacheIterator>
CachingFasta::find_cached(const GenomicRegion& request_region) const noexcept
{
    if (!is_contig_cached(request_region)) return boost::none;
    const auto& contig_cache = sequence_cache_.at(request_region.contig_name());
    if (contig_cache.empty()) return boost::none;
    const auto& contig_region = request_region.contig_region();
    auto iter = contig_cache.lower_bound(contig_region);
    // iter now points to the first region that is not before the request region
    if (iter != std::cbegin(contig_cache)
        && (iter == std::cend(contig_cache) || !begins_equal(contig_region, iter->first))) {
        --iter;
    }
    if (contains(iter->first, contig_region)) {
        return iter;
    } else {
        return boost::none;
    }
}

GenomicRegion CachingFasta::get_region_to_fetch(const GenomicRegion& requested_region) const
{
    // We know the entire region is not in cache, but parts of it may be
    if (sequence_cache_.count(requested_region.contig_name()) == 0) {
        return get_new_contig_chunk(requested_region);
    } else {
        return get_partial_contig_chunk(requested_region);
    }
}

CachingFasta::GenomicSize CachingFasta::get_lhs_extension_size(const GenomicRegion& requested_region) const
{
    assert(max_cache_size_ >= size(requested_region));
    return std::min(requested_region.begin(),
                    static_cast<GenomicSize>((max_cache_size_ - size(requested_region)) * locality_bias_ * (1.0 - forward_bias_)));
}

CachingFasta::GenomicSize CachingFasta::get_rhs_extension_size(const GenomicRegion& requested_region) const
{
    const auto contig_size = contig_sizes_.at(requested_region.contig_name());
    const auto remaining_rhs_size = contig_size - std::min(requested_region.end(), contig_size);
    assert(max_cache_size_ >= size(requested_region));
    const auto remaining_cache_size = max_cache_size_ - size(requested_region);
    assert(locality_bias_ * forward_bias_ <= 1.0);
    const auto preferred_extension_size = static_cast<GenomicSize>(remaining_cache_size * locality_bias_ * forward_bias_);
    return std::min(remaining_rhs_size, preferred_extension_size);
}

GenomicRegion CachingFasta::get_new_contig_chunk(const GenomicRegion& requested_region) const
{
    return expand(requested_region,
                  get_lhs_extension_size(requested_region),
                  get_rhs_extension_size(requested_region));
}

GenomicRegion CachingFasta::get_partial_contig_chunk(const GenomicRegion& requested_region) const
{
    // TODO: optimise this to avoid requesting sequence we already have in cache
    return get_new_contig_chunk(requested_region);
}

void CachingFasta::add_sequence_to_cache(GeneticSequence&& sequence, GenomicRegion&& region) const
{
    assert(size(region) <= max_cache_size_);
    recache_overlapped_regions(sequence, region);
    const auto& contig = region.contig_name();
    if (sequence_cache_.count(contig) == 0 && contig_sizes_.at(contig) >= get_remaining_cache_size()) {
        // Try to clear some room for the new contig hit as we can expect more
        // hits to come
        const auto target_cache_size = static_cast<GenomicSize>(current_cache_size_ * (1.0 - locality_bias_));
        reduce_cache(target_cache_size);
    }
    sequence_cache_[contig].emplace(region.contig_region(), std::move(sequence));
    current_cache_size_ += size(region);
    recently_used_regions_.push_front(std::move(region));
    reduce_cache(max_cache_size_);
    assert(!recently_used_regions_.empty());
}

void CachingFasta::register_cache_hit(const GenomicRegion& region) const
{
    const auto hit = std::find_if(std::cbegin(recently_used_regions_), std::cend(recently_used_regions_),
                                 [&region] (const auto& cached_region) {
                                     return contains(cached_region, region);
                                 });
    assert(hit != std::cend(recently_used_regions_));
    if (hit != std::cbegin(recently_used_regions_)) {
        recently_used_regions_.splice(std::begin(recently_used_regions_),
                                      recently_used_regions_, hit, std::next(hit));
    }
}

CachingFasta::OverlapRange CachingFasta::overlap_range(const GenomicRegion& region) const
{
    assert(sequence_cache_.count(region.contig_name()) == 1);
    const auto& contig_cache = sequence_cache_.at(region.contig_name());
    const auto& contig_region = region.contig_region();
    auto p = contig_cache.equal_range(contig_region);
    // p.first now points to the first region that is not before the request region
    if (p.first != std::cbegin(contig_cache)
        && (p.first == std::cend(contig_cache) || !begins_equal(contig_region, p.first->first))) {
        --p.first;
        if (p.first->first.end() <= contig_region.begin()) {
            ++p.first;
        }
    }
    // p.second points to the first region that is greater than the request region
    if (p.second != std::cend(contig_cache) && contig_region.end() > p.second->first.begin()) {
        ++p.second;
    }
    return {p.first, p.second};
}

void CachingFasta::remove_from_sequence_cache(const GenomicRegion& region) const
{
    const auto& contig = region.contig_name();
    const auto iter = sequence_cache_.find(contig);
    assert(iter != std::cend(sequence_cache_));
    iter->second.erase(region.contig_region());
    if (iter->second.empty()) {
        sequence_cache_.erase(iter);
    }
}

void CachingFasta::remove_from_usage_cache(const GenomicRegion& region) const
{
    recently_used_regions_.remove(region);
}

void CachingFasta::replace_in_usage_cache(const GenomicRegion& from, const GenomicRegion& to) const
{
    const auto it = std::find(std::begin(recently_used_regions_), std::end(recently_used_regions_), from);
    assert(it != std::end(recently_used_regions_));
    *it = to;
}

auto get_nonoverlapped_part(const GenomicRegion& lhs, const GenomicRegion& rhs)
{
    return begins_before(lhs, rhs) ? left_overhang_region(lhs, rhs) : right_overhang_region(lhs, rhs);
}

void CachingFasta::recache_overlapped_regions(const GeneticSequence& sequence, const GenomicRegion& region) const
{
    if (sequence_cache_.count(region.contig_name()) == 1) {
        auto cached_overlapped = this->overlap_range(region);
        const auto num_overlapped = std::distance(std::cbegin(cached_overlapped), std::cend(cached_overlapped));
        
        if (num_overlapped > 0) {
            std::vector<ContigSequenceCache::value_type> to_recache {};
            to_recache.reserve(num_overlapped);
            
            std::for_each(std::begin(cached_overlapped), std::end(cached_overlapped),
                          [this, &region, &to_recache] (const ContigSequenceCache::value_type& p) {
                const ContigRegion& cached_contig_region {p.first};
                const GenomicRegion cached_region {region.contig_name(), cached_contig_region};
                current_cache_size_ -= size(cached_contig_region);
                if (contains(region.contig_region(), cached_contig_region)) {
                    remove_from_usage_cache(cached_region);
                } else {
                    const auto cached_subregion_to_keep = get_nonoverlapped_part(cached_region, region);
                    assert(contains(cached_region, cached_subregion_to_keep));
                    to_recache.emplace_back(cached_subregion_to_keep.contig_region(),
                                            get_subsequence(cached_subregion_to_keep.contig_region(),
                                                            cached_contig_region,
                                                            p.second));
                    replace_in_usage_cache(cached_region, cached_subregion_to_keep);
                    current_cache_size_ += size(cached_subregion_to_keep);
                }
            });
            assert(current_cache_size_ <= max_cache_size_);
            // update sequence_cache_ here as std::map iterators are invalidated on erase
            sequence_cache_.at(region.contig_name()).erase(std::begin(cached_overlapped),
                                                           std::end(cached_overlapped));
            sequence_cache_.at(region.contig_name()).insert(std::make_move_iterator(std::begin(to_recache)),
                                                            std::make_move_iterator(std::end(to_recache)));
        }
    }
}

void CachingFasta::reduce_cache(const std::size_t target_size) const
{
    while (current_cache_size_ > target_size) {
        assert(!recently_used_regions_.empty());
        remove_from_sequence_cache(recently_used_regions_.back());
        current_cache_size_ -= size(recently_used_regions_.back());
        recently_used_regions_.pop_back();
    }
}

CachingFasta::GeneticSequence CachingFasta::get_subsequence(const ContigRegion& requested_region,
                                                            const ContigRegion& sequence_region,
                                                            const GeneticSequence& sequence) const
{
    assert(size(sequence_region) == sequence.size());
    assert(contains(sequence_region, requested_region));
    return sequence.substr(begin_distance(sequence_region, requested_region), size(requested_region));
}

} // namespace io
} // namespace octopus
