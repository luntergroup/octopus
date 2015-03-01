//
//  fasta.cpp
//  Octopus
//
//  Created by Daniel Cooke on 10/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "fasta.h"

#include <stdexcept>

#include "genomic_region.h"

Fasta::Fasta(std::string fasta_path)
:Fasta {fasta_path, fasta_path + ".fai"}
{}

Fasta::Fasta(std::string fasta_path, std::string fasta_index_path)
:fasta_path_ {std::move(fasta_path)},
 fasta_index_path_ {std::move(fasta_index_path)}
{
    if (!fs::exists(fasta_path_)) {
        throw std::runtime_error {"Cannot find FASTA " + fasta_path};
    }
    if (!fs::exists(fasta_index_path_)) {
        fasta_index_path_ = fasta_path_.replace_extension("fai");
        if (!fs::exists(fasta_index_path_)) {
            throw std::runtime_error {"Cannot find FASTA index"};
        }
    }
    if (!is_valid_fasta()) {
        throw std::runtime_error {"Invalid FASTA " + fasta_path};
    }
    
    fasta_ = std::ifstream(fasta_path_.string());
    fasta_contig_indices_ = bioio::read_fasta_index(fasta_index_path_.string());
}

bool Fasta::is_valid_fasta() const
{
    auto extension = fasta_path_.extension().string();
    if (extension != ".fa" && extension != ".fasta") {
        return false;
    }
    return true; // TODO: could actually check valid fasta format
}

std::string Fasta::get_reference_name()
{
    return fasta_path_.stem().string();
}

std::vector<std::string> Fasta::get_contig_names()
{
    return bioio::get_fasta_index_contig_names(fasta_index_path_.string());
}

Fasta::SizeType Fasta::get_contig_size(std::string contig_name)
{
    return static_cast<SizeType>(fasta_contig_indices_.at(contig_name).length);
}

Fasta::SequenceType Fasta::get_sequence(const GenomicRegion& a_region)
{
    return bioio::read_fasta_contig(fasta_, fasta_contig_indices_.at(a_region.get_contig_name()),
                                    a_region.get_begin(), size(a_region));
}
