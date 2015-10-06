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

namespace fs = boost::filesystem;

class GenomicRegion;

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
    ~CachingFasta() noexcept override = default;
    
    CachingFasta(const CachingFasta&)            = default;
    CachingFasta& operator=(const CachingFasta&) = default;
    CachingFasta(CachingFasta&&)                 = default;
    CachingFasta& operator=(CachingFasta&&)      = default;
    
private:
    Fasta fasta_;
    
    using ContigSequenceCache = std::map<GenomicRegion, SequenceType>;
    using CacheIterator       = ContigSequenceCache::const_iterator;
    using OverlapRange        = boost::iterator_range<CacheIterator>;
    
    std::string do_get_reference_name() const override;
    std::vector<std::string> do_get_contig_names() override;
    SizeType do_get_contig_size(const std::string& contig_name) override;
    SequenceType do_fetch_sequence(const GenomicRegion& region) override;
    
    std::unordered_map<std::string, SizeType> contig_size_cache_;
    
    std::unordered_map<std::string, ContigSequenceCache> sequence_cache_;
    std::list<GenomicRegion> recently_used_regions_;
    
    SizeType used_cache_size_      = 0;
    const SizeType max_cache_size_ = 10'000'000;
    
    void setup_cache();
    
    GenomicRegion get_region_to_fetch(const GenomicRegion& requested_region) const;
    
    GenomicRegion get_new_contig_chunk(const GenomicRegion& requested_region) const;
    GenomicRegion get_hit_contig_chunk(const GenomicRegion& requested_region) const;
    
    bool is_region_cached(const GenomicRegion& region) const;
    void add_sequence_to_cache(const SequenceType& sequence, const GenomicRegion& region);
    void update_cache_position(const GenomicRegion& region);
    CacheIterator get_cache_iterator(const GenomicRegion& requested_region) const;
    OverlapRange overlap_range(const GenomicRegion& region) const;
    void remove_from_cache(const GenomicRegion& region);
    void recache_overlapped_regions(const SequenceType& sequence, const GenomicRegion& region);
    
    SequenceType get_subsequence(const GenomicRegion& requested_region,
                                 const GenomicRegion& sequence_region,
                                 const SequenceType& sequence) const;
};

#endif /* defined(__Octopus__caching_fasta__) */
