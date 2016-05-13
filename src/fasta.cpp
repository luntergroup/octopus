//
//  fasta.cpp
//  Octopus
//
//  Created by Daniel Cooke on 10/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "fasta.hpp"

#include <iostream>
#include <stdexcept>

#include <boost/filesystem/operations.hpp>

#include "genomic_region.hpp"

Fasta::Fasta(Path fasta_path)
:
Fasta {fasta_path, fasta_path.string() + ".fai"}
{}

Fasta::Fasta(Path fasta_path, Path fasta_index_path)
:
fasta_path_ {std::move(fasta_path)},
fasta_index_path_ {std::move(fasta_index_path)}
{
    using boost::filesystem::exists;
    
    if (!exists(fasta_path_)) {
        throw std::runtime_error {"Fasta: given fasta path " + fasta_path_.string() + " does not exist"};
    }
    
    if (!exists(fasta_index_path_)) {
        fasta_index_path_ = fasta_path_.replace_extension("fai");
        
        if (!exists(fasta_index_path_)) {
            throw std::runtime_error { "Fasta: given fasta index path " + fasta_index_path_.string()
                + " does not exist"};
        }
    }
    
    if (is_valid()) {
        fasta_ = std::ifstream(fasta_path_.string());
        fasta_index_ = bioio::read_fasta_index(fasta_index_path_.string());
    }
}

// virtual private methods

bool Fasta::do_is_open() const noexcept
{
    try {
        return fasta_.is_open();
    } catch (...) {
        // only because std::ifstream::is_open is not declared noexcept
        return false;
    }
}

std::string Fasta::do_fetch_reference_name() const
{
    return fasta_path_.stem().string();
}

std::vector<Fasta::ContigNameType> Fasta::do_fetch_contig_names() const
{
    return bioio::read_fasta_index_contig_names(fasta_index_path_.string());
}

Fasta::SizeType Fasta::do_fetch_contig_size(const ContigNameType& contig) const
{
    if (fasta_index_.count(contig) == 0) {
        throw std::runtime_error {"contig \"" + contig +
            "\" not found in fasta index \"" + fasta_index_path_.string() + "\""};
    }
    return static_cast<SizeType>(fasta_index_.at(contig).length);
}

Fasta::SequenceType Fasta::do_fetch_sequence(const GenomicRegion& region) const
{
    return bioio::read_fasta_contig(fasta_, fasta_index_.at(contig_name(region)),
                                    mapped_begin(region), size(region));
}

bool Fasta::is_valid() const noexcept
{
    const auto extension = fasta_path_.extension().string();
    
    if (extension != ".fa" && extension != ".fasta") {
        return false;
    }
    
    return true; // TODO: could actually check valid fasta format
}
