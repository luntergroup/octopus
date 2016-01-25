//
//  threadsafe_caching_fasta.hpp
//  Octopus
//
//  Created by Daniel Cooke on 07/09/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__threadsafe_caching_fasta__
#define __Octopus__threadsafe_caching_fasta__

#include <string>
#include <vector>
#include <mutex>

#include <boost/filesystem/path.hpp>

#include "reference_genome_impl.hpp"
#include "caching_fasta.hpp"

class GenomicRegion;

class ThreadsafeCachingFasta : public ReferenceGenomeImpl
{
public:
    using Path = Fasta::Path;
    
    using ContigNameType = ReferenceGenomeImpl::ContigNameType;
    using SizeType       = ReferenceGenomeImpl::SizeType;
    using SequenceType   = ReferenceGenomeImpl::SequenceType;
    
    ThreadsafeCachingFasta() = delete;
    explicit ThreadsafeCachingFasta(Path fasta_path);
    explicit ThreadsafeCachingFasta(Path fasta_path, SizeType max_cache_size);
    explicit ThreadsafeCachingFasta(Path fasta_path, Path fasta_index_path);
    explicit ThreadsafeCachingFasta(Path fasta_path, Path fasta_index_path, SizeType max_cache_size);
    explicit ThreadsafeCachingFasta(Path fasta_path, SizeType max_cache_size, double locality_bias, double forward_bias);
    ~ThreadsafeCachingFasta() noexcept override = default;
    
    ThreadsafeCachingFasta(const ThreadsafeCachingFasta&)            = default;
    ThreadsafeCachingFasta& operator=(const ThreadsafeCachingFasta&) = default;
    ThreadsafeCachingFasta(ThreadsafeCachingFasta&&)                 = default;
    ThreadsafeCachingFasta& operator=(ThreadsafeCachingFasta&&)      = default;
    
private:
    CachingFasta fasta_;
    mutable std::mutex fasta_mutex_;
    
    bool do_is_open() const noexcept override;
    std::string do_get_reference_name() const override;
    std::vector<ContigNameType> do_get_contig_names() const override;
    SizeType do_get_contig_size(const ContigNameType& contig) const override;
    SequenceType do_fetch_sequence(const GenomicRegion& region) const override;
};

#endif /* defined(__Octopus__threadsafe_caching_fasta__) */
