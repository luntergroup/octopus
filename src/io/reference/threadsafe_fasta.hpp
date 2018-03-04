// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef threadsafe_fasta_hpp
#define threadsafe_fasta_hpp

#include <string>
#include <vector>
#include <mutex>
#include <memory>

#include <boost/filesystem/path.hpp>

#include "reference_reader.hpp"
#include "fasta.hpp"

namespace octopus {
    
class GenomicRegion;

namespace io {

class ThreadsafeFasta : public ReferenceReader
{
public:
    using Path = Fasta::Path;
    
    using ContigName       = ReferenceReader::ContigName;
    using GenomicSize      = ReferenceReader::GenomicSize;
    using GeneticSequence  = ReferenceReader::GeneticSequence;
    
    ThreadsafeFasta() = delete;
    
    ThreadsafeFasta(std::unique_ptr<Fasta> fasta);
    
    ThreadsafeFasta(const ThreadsafeFasta&);
    ThreadsafeFasta& operator=(ThreadsafeFasta);
    ThreadsafeFasta(ThreadsafeFasta&&);
    ThreadsafeFasta& operator=(ThreadsafeFasta&&);
    
private:
    std::unique_ptr<Fasta> fasta_;
    mutable std::mutex mutex_;
    
    std::unique_ptr<ReferenceReader> do_clone() const override;
    bool do_is_open() const noexcept override;
    std::string do_fetch_reference_name() const override;
    std::vector<ContigName> do_fetch_contig_names() const override;
    GenomicSize do_fetch_contig_size(const ContigName& contig) const override;
    GeneticSequence do_fetch_sequence(const GenomicRegion& region) const override;
};
    
} // namespace io
} // namespace octopus

#endif
