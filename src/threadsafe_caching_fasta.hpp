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
    
    using SequenceType = ReferenceGenomeImpl::SequenceType;
    using SizeType     = ReferenceGenomeImpl::SizeType;
    
    ThreadsafeCachingFasta() = delete;
    explicit ThreadsafeCachingFasta(Path fasta_path);
    explicit ThreadsafeCachingFasta(Path fasta_path, Path fasta_index_path);
    ~ThreadsafeCachingFasta() noexcept override = default;
    
    ThreadsafeCachingFasta(const ThreadsafeCachingFasta&)            = default;
    ThreadsafeCachingFasta& operator=(const ThreadsafeCachingFasta&) = default;
    ThreadsafeCachingFasta(ThreadsafeCachingFasta&&)                 = default;
    ThreadsafeCachingFasta& operator=(ThreadsafeCachingFasta&&)      = default;
    
private:
    CachingFasta fasta_;
    mutable std::mutex fasta_mutex_;
    
    std::string do_get_reference_name() const override;
    std::vector<std::string> do_get_contig_names() const override;
    SizeType do_get_contig_size(const std::string& contig_name) const override;
    SequenceType do_fetch_sequence(const GenomicRegion& region) const override;
};

#endif /* defined(__Octopus__threadsafe_caching_fasta__) */
