// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "threadsafe_fasta.hpp"

#include <basics/genomic_region.hpp>

namespace octopus { namespace io {

ThreadsafeFasta::ThreadsafeFasta(std::unique_ptr<Fasta> fasta)
: fasta_ {std::move(fasta)}
{}

ThreadsafeFasta::ThreadsafeFasta(const ThreadsafeFasta& other)
{
    std::lock_guard<std::mutex> lock {other.mutex_};
    fasta_ = std::make_unique<Fasta>(*other.fasta_);
}

ThreadsafeFasta& ThreadsafeFasta::operator=(ThreadsafeFasta other)
{
    std::swap(fasta_, other.fasta_);
    return *this;
}

ThreadsafeFasta::ThreadsafeFasta(ThreadsafeFasta&& other)
{
    std::lock_guard<std::mutex> lock {other.mutex_};
    fasta_ = std::move(other.fasta_);
}

std::unique_ptr<ReferenceReader> ThreadsafeFasta::do_clone() const
{
    return std::make_unique<ThreadsafeFasta>(*this);
}

bool ThreadsafeFasta::do_is_open() const noexcept
{
    return fasta_->is_open();
}

std::string ThreadsafeFasta::do_fetch_reference_name() const
{
    return fasta_->fetch_reference_name(); // don't need mutex as const
}

std::vector<ThreadsafeFasta::ContigName> ThreadsafeFasta::do_fetch_contig_names() const
{
    std::lock_guard<std::mutex> lock {mutex_};
    return fasta_->fetch_contig_names();
}

ThreadsafeFasta::GenomicSize ThreadsafeFasta::do_fetch_contig_size(const ContigName& contig) const
{
    std::lock_guard<std::mutex> lock {mutex_};
    return fasta_->fetch_contig_size(contig);
}

ThreadsafeFasta::GeneticSequence ThreadsafeFasta::do_fetch_sequence(const GenomicRegion& region) const
{
    std::lock_guard<std::mutex> lock {mutex_};
    return fasta_->fetch_sequence(region);
}

} // namespace io
} // namespace octopus
