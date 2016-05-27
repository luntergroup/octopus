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

class ThreadsafeFasta : public ReferenceGenomeImpl
{
public:
    using Path = Fasta::Path;
    
    using ContigNameType = ReferenceGenomeImpl::ContigNameType;
    using SizeType       = ReferenceGenomeImpl::SizeType;
    using SequenceType   = ReferenceGenomeImpl::SequenceType;
    
    ThreadsafeFasta() = delete;
    
    ThreadsafeFasta(std::unique_ptr<Fasta> fasta);
    
    ~ThreadsafeFasta() noexcept override = default;
    
    ThreadsafeFasta(const ThreadsafeFasta&)            = default;
    ThreadsafeFasta& operator=(const ThreadsafeFasta&) = default;
    ThreadsafeFasta(ThreadsafeFasta&&)                 = default;
    ThreadsafeFasta& operator=(ThreadsafeFasta&&)      = default;
    
private:
    std::unique_ptr<Fasta> fasta_;
    mutable std::mutex mutex_;
    
    bool do_is_open() const noexcept override;
    std::string do_fetch_reference_name() const override;
    std::vector<ContigNameType> do_fetch_contig_names() const override;
    SizeType do_fetch_contig_size(const ContigNameType& contig) const override;
    SequenceType do_fetch_sequence(const GenomicRegion& region) const override;
};

#endif /* defined(__Octopus__threadsafe_fasta__) */
