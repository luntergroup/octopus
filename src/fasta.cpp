//
//  fasta.cpp
//  Octopus
//
//  Created by Daniel Cooke on 10/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "fasta.hpp"

#include <stdexcept>

#include <boost/filesystem/operations.hpp>

#include "genomic_region.hpp"

Fasta::Fasta(fs::path fasta_path)
:
Fasta {fasta_path, fasta_path.string() + ".fai"}
{}

Fasta::Fasta(fs::path fasta_path, fs::path fasta_index_path)
:
fasta_path_ {std::move(fasta_path)},
fasta_index_path_ {std::move(fasta_index_path)}
{
    if (!fs::exists(fasta_path_)) {
        throw std::runtime_error {"Cannot find FASTA \"" + fasta_path.string() + "\""};
    }
    
    if (!fs::exists(fasta_index_path_)) {
        fasta_index_path_ = fasta_path_.replace_extension("fai");
        
        if (!fs::exists(fasta_index_path_)) {
            throw std::runtime_error {"Cannot find FASTA index \"" + fasta_index_path_.string() + "\""};
        }
    }
    
    if (!is_valid_fasta()) {
        throw std::runtime_error {"Invalid FASTA \"" + fasta_path.string() + "\""};
    }
    
    fasta_ = std::ifstream(fasta_path_.string());
    fasta_contig_indices_ = bioio::read_fasta_index(fasta_index_path_.string());
}

std::string Fasta::do_get_reference_name() const
{
    return fasta_path_.stem().string();
}

std::vector<std::string> Fasta::do_get_contig_names() const
{
    return bioio::get_fasta_index_contig_names(fasta_index_path_.string());
}

Fasta::SizeType Fasta::do_get_contig_size(const std::string& contig_name) const
{
    if (fasta_contig_indices_.count(contig_name) == 0) {
        throw std::runtime_error {"contig \"" + contig_name +
            "\" not found in fasta index \"" + fasta_index_path_.string() + "\""};
    }
    
    return static_cast<SizeType>(fasta_contig_indices_.at(contig_name).length);
}

Fasta::SequenceType Fasta::do_fetch_sequence(const GenomicRegion& region) const
{
    return bioio::read_fasta_contig(fasta_, fasta_contig_indices_.at(get_contig_name(region)),
                                    get_begin(region), size(region));
}

bool Fasta::is_valid_fasta() const noexcept
{
    const auto extension = fasta_path_.extension().string();
    
    if (extension != ".fa" && extension != ".fasta") {
        return false;
    }
    
    return true; // TODO: could actually check valid fasta format
}
