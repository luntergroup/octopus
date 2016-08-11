// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef caching_fasta_hpp
#define caching_fasta_hpp

#include <unordered_map>
#include <map>
#include <list>
#include <cstddef>
#include <mutex>
#include <memory>

#include <boost/filesystem/path.hpp>
#include <boost/range/iterator_range_core.hpp>

#include <basics/contig_region.hpp>

#include "reference_reader.hpp"
#include "fasta.hpp"

namespace octopus {
    
class GenomicRegion;

namespace io {

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

class CachingFasta : public ReferenceReader
{
public:
    using Path = Fasta::Path;
    
    using ContigName      = ReferenceReader::ContigName;
    using GenomicSize     = ReferenceReader::GenomicSize;
    using GeneticSequence = ReferenceReader::GeneticSequence;
    
    CachingFasta() = delete;
    
    CachingFasta(std::unique_ptr<ReferenceReader> fasta);
    CachingFasta(std::unique_ptr<ReferenceReader> fasta, GenomicSize max_cache_size);
    CachingFasta(std::unique_ptr<ReferenceReader> fasta, GenomicSize max_cache_size,
                 double locality_bias, double forward_bias);
    
    CachingFasta(const CachingFasta&);
    CachingFasta& operator=(CachingFasta);
    CachingFasta(CachingFasta&&);
    CachingFasta& operator=(CachingFasta&&)      = default;
    
private:
    using ContigSequenceCache = std::map<ContigRegion, GeneticSequence>;
    using SequenceCache       = std::unordered_map<std::string, ContigSequenceCache>;
    using CacheIterator       = ContigSequenceCache::const_iterator;
    using OverlapRange        = boost::iterator_range<CacheIterator>;
    
    std::unique_ptr<ReferenceReader> fasta_;
    
    std::unordered_map<ContigName, GenomicSize> contig_sizes_;
    
    mutable SequenceCache sequence_cache_;
    mutable std::list<GenomicRegion> recently_used_regions_; // TODO: also make into map of contigs?
    
    GenomicSize genome_size_;
    
    GenomicSize max_cache_size_;
    mutable GenomicSize current_cache_size_;
    
    double locality_bias_;
    double forward_bias_;
    
    mutable std::mutex mutex_;
    
    std::unique_ptr<ReferenceReader> do_clone() const override;
    bool do_is_open() const noexcept override;
    std::string do_fetch_reference_name() const override;
    std::vector<ContigName> do_fetch_contig_names() const override;
    GenomicSize do_fetch_contig_size(const ContigName& contig) const override;
    GeneticSequence do_fetch_sequence(const GenomicRegion& region) const override;
    
    void setup_cache();
    GenomicSize get_remaining_cache_size() const;
    GenomicSize get_lhs_extension_size(const GenomicRegion& requested_region) const;
    GenomicSize get_rhs_extension_size(const GenomicRegion& requested_region) const;
    GenomicRegion get_region_to_fetch(const GenomicRegion& requested_region) const;
    GenomicRegion get_new_contig_chunk(const GenomicRegion& requested_region) const;
    GenomicRegion get_partial_contig_chunk(const GenomicRegion& requested_region) const;
    bool is_contig_cached(const GenomicRegion& region) const;
    CacheIterator find_cache_iterator(const GenomicRegion& requested_region) const;
    void add_sequence_to_cache(GeneticSequence&& sequence, GenomicRegion&& region) const;
    void update_cache_position(const GenomicRegion& region) const;
    OverlapRange overlap_range(const GenomicRegion& region) const;
    void remove_from_sequence_cache(const GenomicRegion& region) const;
    void remove_from_usage_cache(const GenomicRegion& region) const;
    void replace_in_usage_cache(const GenomicRegion& from, const GenomicRegion& to) const;
    void recache_overlapped_regions(const GeneticSequence& sequence, const GenomicRegion& region) const;
    GeneticSequence get_subsequence(const ContigRegion& requested_region,
                                    const ContigRegion& sequence_region,
                                    const GeneticSequence& sequence) const;
};

} // namespace io
} // namespace octopus

#endif
