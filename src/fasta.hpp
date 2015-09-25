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
#include <unordered_map>
#include <boost/filesystem/path.hpp>

#include "reference_genome_impl.hpp"
#include "bioio.hpp"

class GenomicRegion;

namespace fs = boost::filesystem;

class Fasta : public ReferenceGenomeImpl
{
public:
    using SequenceType = ReferenceGenomeImpl::SequenceType;
    using SizeType     = ReferenceGenomeImpl::SizeType;
    
    Fasta() = delete;
    explicit Fasta(fs::path fasta_path);
    explicit Fasta(fs::path fasta_path, fs::path fasta_index_path);
    ~Fasta() noexcept override = default;
    
    Fasta(const Fasta&)            = default;
    Fasta& operator=(const Fasta&) = default;
    Fasta(Fasta&&)                 = default;
    Fasta& operator=(Fasta&&)      = default;

private:
    fs::path fasta_path_;
    fs::path fasta_index_path_;
    std::ifstream fasta_;
    std::unordered_map<std::string, bioio::FastaIndex> fasta_contig_indices_;
    
    std::string do_get_reference_name() const override;
    std::vector<std::string> do_get_contig_names() override;
    SizeType do_get_contig_size(const std::string& contig_name) override;
    SequenceType do_get_sequence(const GenomicRegion& region) override;
    
    bool is_valid_fasta() const noexcept;
};

#endif /* defined(__Octopus__fasta__) */
