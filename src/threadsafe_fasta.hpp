//
//  threadsafe_fasta.hpp
//  Octopus
//
//  Created by Daniel Cooke on 17/08/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__threadsafe_fasta__
#define __Octopus__threadsafe_fasta__

#include <string>
#include <vector>
#include <mutex>
#include <boost/filesystem/path.hpp>

#include "reference_genome_impl.hpp"
#include "fasta.hpp"

class GenomicRegion;

namespace fs = boost::filesystem;

class ThreadsafeFasta : public ReferenceGenomeImpl
{
public:
    using SequenceType = ReferenceGenomeImpl::SequenceType;
    using SizeType     = ReferenceGenomeImpl::SizeType;
    
    ThreadsafeFasta() = delete;
    explicit ThreadsafeFasta(fs::path fasta_path);
    explicit ThreadsafeFasta(fs::path fasta_path, fs::path fasta_index_path);
    ~ThreadsafeFasta() noexcept override = default;
    
    ThreadsafeFasta(const ThreadsafeFasta&)            = default;
    ThreadsafeFasta& operator=(const ThreadsafeFasta&) = default;
    ThreadsafeFasta(ThreadsafeFasta&&)                 = default;
    ThreadsafeFasta& operator=(ThreadsafeFasta&&)      = default;
    
private:
    Fasta fasta_;
    std::mutex mutex_;
    
    std::string do_get_reference_name() const override;
    std::vector<std::string> do_get_contig_names() override;
    SizeType do_get_contig_size(const std::string& contig_name) override;
    SequenceType do_fetch_sequence(const GenomicRegion& region) override;
};

#endif /* defined(__Octopus__threadsafe_fasta__) */
