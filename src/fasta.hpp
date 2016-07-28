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
#include <stdexcept>

#include <boost/filesystem/path.hpp>

#include "reference_genome_impl.hpp"
#include "bioio.hpp"

class GenomicRegion;

class Fasta : public ReferenceGenomeImpl
{
public:
    class MissingFastaIndex;
    
    using Path = boost::filesystem::path;
    
    using ContigName      = ReferenceGenomeImpl::ContigName;
    using GenomicSize     = ReferenceGenomeImpl::GenomicSize;
    using GeneticSequence = ReferenceGenomeImpl::GeneticSequence;
    
    Fasta() = delete;
    
    Fasta(Path fasta_path);
    Fasta(Path fasta_path, Path fasta_index_path);
    
    Fasta(const Fasta&)            = default;
    Fasta& operator=(const Fasta&) = default;
    Fasta(Fasta&&)                 = default;
    Fasta& operator=(Fasta&&)      = default;
    
    ~Fasta() noexcept override = default;

private:
    Path fasta_path_;
    Path fasta_index_path_;
    
    mutable std::ifstream fasta_;
    bioio::FastaIndex fasta_index_;
    
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

#endif /* defined(__Octopus__fasta__) */
