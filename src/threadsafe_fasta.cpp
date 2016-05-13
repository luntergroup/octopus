//
//  threadsafe_fasta.cpp
//  Octopus
//
//  Created by Daniel Cooke on 17/08/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "threadsafe_fasta.hpp"

#include "genomic_region.hpp"

ThreadsafeFasta::ThreadsafeFasta(std::unique_ptr<Fasta> fasta)
:
fasta_ {std::move(fasta)}
{}

bool ThreadsafeFasta::do_is_open() const noexcept
{
    return fasta_->is_open();
}

std::string ThreadsafeFasta::do_fetch_reference_name() const
{
    return fasta_->fetch_reference_name(); // don't need mutex as const
}

std::vector<ThreadsafeFasta::ContigNameType> ThreadsafeFasta::do_fetch_contig_names() const
{
    std::lock_guard<std::mutex> lock {mutex_};
    return fasta_->fetch_contig_names();
}

ThreadsafeFasta::SizeType ThreadsafeFasta::do_fetch_contig_size(const ContigNameType& contig) const
{
    std::lock_guard<std::mutex> lock {mutex_};
    return fasta_->fetch_contig_size(contig);
}

ThreadsafeFasta::SequenceType ThreadsafeFasta::do_fetch_sequence(const GenomicRegion& region) const
{
    std::lock_guard<std::mutex> lock {mutex_};
    return fasta_->fetch_sequence(region);
}
