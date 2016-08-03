//
//  fasta.hpp
//  Octopus
//
//  Created by Daniel Cooke on 10/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__fasta__
#define __Octopus__fasta__

#include <string>
#include <vector>
#include <cstdint>
#include <fstream>
#include <memory>
#include <stdexcept>

#include <boost/filesystem/path.hpp>

#include <bioio.hpp>

#include "reference_reader.hpp"

namespace octopus {
    
class GenomicRegion;

namespace io {

class Fasta : public ReferenceReader
{
public:
    class MissingFastaIndex;
    
    using Path = boost::filesystem::path;
    
    using ContigName      = ReferenceReader::ContigName;
    using GenomicSize     = ReferenceReader::GenomicSize;
    using GeneticSequence = ReferenceReader::GeneticSequence;
    
    Fasta() = delete;
    
    Fasta(Path fasta_path);
    Fasta(Path fasta_path, Path fasta_index_path);
    
    Fasta(const Fasta&);
    Fasta& operator=(Fasta);
    Fasta(Fasta&&)            = default;
    Fasta& operator=(Fasta&&) = default;
    
private:
    Path path_;
    Path index_path_;
    
    mutable std::ifstream fasta_;
    bioio::FastaIndex fasta_index_;
    
    std::unique_ptr<ReferenceReader> do_clone() const override;
    bool do_is_open() const noexcept override;
    std::string do_fetch_reference_name() const override;
    std::vector<ContigName> do_fetch_contig_names() const override;
    GenomicSize do_fetch_contig_size(const ContigName& contig) const override;
    GeneticSequence do_fetch_sequence(const GenomicRegion& region) const override;
    
    bool is_valid() const noexcept;
};

class Fasta::MissingFastaIndex : public std::runtime_error
{
public:
    using Path = Fasta::Path;
    
    MissingFastaIndex(Path fasta_path, Path expected_index_path);
    
    virtual ~MissingFastaIndex() noexcept = default;
    
    const Path& fasta_path() const noexcept;
    const Path& expected_index_path() const noexcept;
    
    virtual const char* what() const noexcept;
    
private:
    Path fasta_path_, expected_index_path_;
    mutable std::string msg_;
};

} // namespace io
} // namespace octopus

#endif /* defined(__Octopus__fasta__) */
