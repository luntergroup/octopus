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
#include <mutex>

#include <boost/filesystem/path.hpp>
#include <boost/range/iterator_range_core.hpp>

#include "reference_genome_impl.hpp"
#include "fasta.hpp"
#include "contig_region.hpp"
#include "genomic_region.hpp"

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
    using Path = Fasta::Path;
    
    using ContigNameType = ReferenceGenomeImpl::ContigNameType;
    using SizeType       = ReferenceGenomeImpl::SizeType;
    using SequenceType   = ReferenceGenomeImpl::SequenceType;
    
    CachingFasta() = delete;
    
    explicit CachingFasta(std::unique_ptr<ReferenceGenomeImpl> fasta);
    explicit CachingFasta(std::unique_ptr<ReferenceGenomeImpl> fasta, SizeType max_cache_size);
    explicit CachingFasta(std::unique_ptr<ReferenceGenomeImpl> fasta,  SizeType max_cache_size,
                          double locality_bias, double forward_bias);
    
    ~CachingFasta() noexcept override = default;
    
    CachingFasta(const CachingFasta&)            = default;
    CachingFasta& operator=(const CachingFasta&) = default;
    CachingFasta(CachingFasta&&)                 = default;
    CachingFasta& operator=(CachingFasta&&)      = default;
    
private:
    std::unique_ptr<ReferenceGenomeImpl> fasta_;
    
    using ContigSequenceCache = std::map<ContigRegion, SequenceType>;
    using SequenceCache       = std::unordered_map<std::string, ContigSequenceCache>;
    using CacheIterator       = ContigSequenceCache::const_iterator;
    using OverlapRange        = boost::iterator_range<CacheIterator>;
    
    bool do_is_open() const noexcept override;
    std::string do_fetch_reference_name() const override;
    std::vector<ContigNameType> do_fetch_contig_names() const override;
    SizeType do_fetch_contig_size(const ContigNameType& contig) const override;
    SequenceType do_fetch_sequence(const GenomicRegion& region) const override;
    
    std::unordered_map<ContigNameType, SizeType> contig_sizes_;
    
    mutable SequenceCache sequence_cache_;
    mutable std::list<GenomicRegion> recently_used_regions_; // TODO: also make into map of contigs?
    
    SizeType genome_size_;
    
    SizeType max_cache_size_;
    mutable SizeType current_cache_size_;
    
    double locality_bias_;
    double forward_bias_;
    
    mutable std::mutex mutex_;
    
    void setup_cache();
    SizeType get_remaining_cache_size() const;
    SizeType get_lhs_extension_size(const GenomicRegion& requested_region) const;
    SizeType get_rhs_extension_size(const GenomicRegion& requested_region) const;
    GenomicRegion get_region_to_fetch(const GenomicRegion& requested_region) const;
    GenomicRegion get_new_contig_chunk(const GenomicRegion& requested_region) const;
    GenomicRegion get_partial_contig_chunk(const GenomicRegion& requested_region) const;
    bool is_contig_cached(const GenomicRegion& region) const;
    CacheIterator find_cache_iterator(const GenomicRegion& requested_region) const;
    void add_sequence_to_cache(SequenceType&& sequence, GenomicRegion&& region) const;
    void update_cache_position(const GenomicRegion& region) const;
    OverlapRange overlap_range(const GenomicRegion& region) const;
    void remove_from_sequence_cache(const GenomicRegion& region) const;
    void remove_from_usage_cache(const GenomicRegion& region) const;
    void replace_in_usage_cache(const GenomicRegion& from, const GenomicRegion& to) const;
    void recache_overlapped_regions(const SequenceType& sequence, const GenomicRegion& region) const;
    SequenceType get_subsequence(const ContigRegion& requested_region, const ContigRegion& sequence_region,
                                 const SequenceType& sequence) const;
};

#endif /* defined(__Octopus__caching_fasta__) */
