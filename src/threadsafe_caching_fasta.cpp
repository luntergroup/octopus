//
//  threadsafe_caching_fasta.cpp
//  Octopus
//
//  Created by Daniel Cooke on 07/09/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "threadsafe_caching_fasta.hpp"

#include "genomic_region.hpp"

ThreadsafeCachingFasta::ThreadsafeCachingFasta(Path fasta_path)
:
fasta_ {std::move(fasta_path)},
fasta_mutex_ {}
{}

ThreadsafeCachingFasta::ThreadsafeCachingFasta(Path fasta_path, Path fasta_index_path)
:
fasta_ {std::move(fasta_path), std::move(fasta_index_path)},
fasta_mutex_ {}
{}

std::string ThreadsafeCachingFasta::do_get_reference_name() const
{
    return fasta_.get_reference_name(); // don't need mutex as const
}

std::vector<ThreadsafeCachingFasta::ContigNameType> ThreadsafeCachingFasta::do_get_contig_names() const
{
    std::lock_guard<std::mutex> lock {fasta_mutex_};
    return fasta_.get_contig_names();
}

ThreadsafeCachingFasta::SizeType ThreadsafeCachingFasta::do_get_contig_size(const ContigNameType& contig) const
{
    std::lock_guard<std::mutex> lock {fasta_mutex_};
    return fasta_.get_contig_size(contig);
}

ThreadsafeCachingFasta::SequenceType ThreadsafeCachingFasta::do_fetch_sequence(const GenomicRegion& region) const
{
    std::lock_guard<std::mutex> lock {fasta_mutex_};
    return fasta_.fetch_sequence(region);
}
