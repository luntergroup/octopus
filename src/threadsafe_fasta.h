//
//  threadsafe_fasta.h
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

#include "i_reference_genome_impl.h"
#include "fasta.h"

class GenomicRegion;

namespace fs = boost::filesystem;

class ThreadsafeFasta : public IReferenceGenomeImpl
{
public:
    using SequenceType = IReferenceGenomeImpl::SequenceType;
    using SizeType     = IReferenceGenomeImpl::SizeType;
    
    ThreadsafeFasta() = delete;
    explicit ThreadsafeFasta(fs::path fasta_path);
    explicit ThreadsafeFasta(fs::path fasta_path, fs::path fasta_index_path);
    ~ThreadsafeFasta() noexcept override = default;
    
    ThreadsafeFasta(const ThreadsafeFasta&)            = default;
    ThreadsafeFasta& operator=(const ThreadsafeFasta&) = default;
    ThreadsafeFasta(ThreadsafeFasta&&)                 = default;
    ThreadsafeFasta& operator=(ThreadsafeFasta&&)      = default;
    
    std::string get_reference_name() const override;
    std::vector<std::string> get_contig_names() override;
    SizeType get_contig_size(const std::string& contig_name) override;
    SequenceType get_sequence(const GenomicRegion& region) override;
    
private:
    Fasta fasta_;
    std::mutex fasta_mutex_;
};

#endif /* defined(__Octopus__threadsafe_fasta__) */
