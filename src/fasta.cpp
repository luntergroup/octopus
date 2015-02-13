//
//  fasta.cpp
//  Octopus
//
//  Created by Daniel Cooke on 10/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "fasta.h"

Fasta::Fasta(std::string fasta_path)
:   Fasta(fasta_path, fasta_path + ".fai")
{}

Fasta::Fasta(std::string fasta_path, std::string fasta_index_path)
:   fasta_path_ {std::move(fasta_path)},
    fasta_index_path_ {std::move(fasta_index_path)}
{
    if (!is_valid_fasta()) {
        throw std::runtime_error {"Invalid FASTA " + fasta_path};
    }
    
    fasta_ = std::ifstream(fasta_path_);
    fasta_contig_indices_ = std::move(bioio::read_fasta_index(fasta_index_path_));
}

bool Fasta::is_valid_fasta() const
{
//    auto extension = fasta_path_.extension().string();
//    if (!(extension == "fa" || extension == "fasta")) {
//        return false;
//    }
    return true; // TODO: could actually check valid fasta format
}

std::string Fasta::get_reference_name()
{
    //return fasta_path_.filename().string();
    // TODO: get boost to work
    auto name_begin = fasta_path_.find_last_of('/') + 1;
    auto name_size  = fasta_path_.find_last_of('.') - name_begin;
    return fasta_path_.substr(name_begin, name_size);
}

std::vector<std::string> Fasta::get_contig_names()
{
    return bioio::get_fasta_index_contig_names(fasta_index_path_);
}

std::uint_fast32_t Fasta::get_contig_size(std::string contig_name)
{
    return static_cast<uint_fast32_t>(fasta_contig_indices_.at(contig_name).length);
}

std::string Fasta::get_sequence(const GenomicRegion& a_region)
{
    if (is_in_cache(a_region)) {
        return region_cache_.at(a_region);
    } else {
        return bioio::read_fasta_contig(fasta_, fasta_contig_indices_.at(a_region.get_contig_name()),
                                        a_region.get_begin(), size(a_region));
    }
}
