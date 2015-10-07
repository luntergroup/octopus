//
//  caching_fasta.hpp
//  Octopus
//
//  Created by Daniel Cooke on 22/07/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__caching_fasta__
#define __Octopus__caching_fasta__

#include <unordered_map>
#include <map>
#include <list>
#include <cstddef>
#include <boost/filesystem/path.hpp>
#include <boost/range/iterator_range_core.hpp>

#include "reference_genome_impl.hpp"
#include "fasta.hpp"
#include "contig_region.hpp"
#include "genomic_region.hpp"

namespace fs = boost::filesystem;

/*
 CachingFasta attempts to reduce the number of file reads in two ways:
 
    1) Caching recently retreived sequence, up to given maximum size.
 
    2) Reading more sequence than requested. The assumption here is that sequence requests are
       clustered in nearby regions.
 
       There are two parameters that control how these 'extended' regions are predicted
 
       - locality bias: the probability the next request region will be nearby the current request region. 
                        If this value is low then there will be little performance gain from using CachingFasta.
 
       - forward bias: the probability the next request region will be to the right hand side of
                       the current request region.
 */

class CachingFasta : public ReferenceGenomeImpl
{
public:
    using SequenceType = ReferenceGenomeImpl::SequenceType;
    using SizeType     = ReferenceGenomeImpl::SizeType;
    
    CachingFasta() = delete;
    explicit CachingFasta(fs::path fasta_path);
    explicit CachingFasta(fs::path fasta_path, SizeType max_cache_size);
    explicit CachingFasta(fs::path fasta_path, fs::path fasta_index_path);
    explicit CachingFasta(fs::path fasta_path, fs::path fasta_index_path, SizeType max_cache_size);
    explicit CachingFasta(fs::path fasta_path, SizeType max_cache_size, double locality_bias, double forward_bias);
    ~CachingFasta() noexcept override = default;
    
    CachingFasta(const CachingFasta&)            = default;
    CachingFasta& operator=(const CachingFasta&) = default;
    CachingFasta(CachingFasta&&)                 = default;
    CachingFasta& operator=(CachingFasta&&)      = default;
    
private:
    Fasta fasta_;
    
    using ContigSequenceCache = std::map<ContigRegion, SequenceType>;
    using SequenceCache       = std::unordered_map<std::string, ContigSequenceCache>;
    using CacheIterator       = ContigSequenceCache::const_iterator;
    using OverlapRange        = boost::iterator_range<CacheIterator>;
    
    std::string do_get_reference_name() const override;
    std::vector<std::string> do_get_contig_names() override;
    SizeType do_get_contig_size(const std::string& contig_name) override;
    SequenceType do_fetch_sequence(const GenomicRegion& region) override;
    
    std::unordered_map<std::string, SizeType> contig_size_cache_;
    
    SequenceCache sequence_cache_;
    std::list<GenomicRegion> recently_used_regions_; // TODO: also make into map of contigs?
    
    SizeType genome_size_;
    
    const SizeType max_cache_size_ = 10'000'000;
    SizeType current_cache_size_   = 0;
    
    const double locality_bias_ = 0.8;
    const double forward_bias_  = 0.8;
    
    void setup_cache();
    SizeType get_remaining_cache_size() const;
    SizeType get_lhs_extension_size(const GenomicRegion& requested_region) const;
    SizeType get_rhs_extension_size(const GenomicRegion& requested_region) const;
    GenomicRegion get_region_to_fetch(const GenomicRegion& requested_region) const;
    GenomicRegion get_new_contig_chunk(const GenomicRegion& requested_region) const;
    GenomicRegion get_hit_contig_chunk(const GenomicRegion& requested_region) const;
    bool is_region_cached(const GenomicRegion& region) const;
    CacheIterator get_cache_iterator(const GenomicRegion& requested_region) const;
    void add_sequence_to_cache(const SequenceType& sequence, const GenomicRegion& region);
    void update_cache_position(const GenomicRegion& region);
    OverlapRange overlap_range(const GenomicRegion& region) const;
    void remove_from_sequence_cache(const GenomicRegion& region);
    void remove_from_usage_cache(const GenomicRegion& region);
    void replace_in_usage_cache(const GenomicRegion& from, const GenomicRegion& to);
    void recache_overlapped_regions(const SequenceType& sequence, const GenomicRegion& region);
    SequenceType get_subsequence(const ContigRegion& requested_region, const ContigRegion& sequence_region,
                                 const SequenceType& sequence) const;
};

#endif /* defined(__Octopus__caching_fasta__) */
