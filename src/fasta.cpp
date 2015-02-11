//
//  fasta.cpp
//  Octopus
//
//  Created by Daniel Cooke on 10/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "fasta.h"

#include "bioio.h"
#include "utils.h"

Fasta::Fasta(std::string fasta_path)
: fasta_path_ {fasta_path}, fasta_index_path_ {}
{
    if (!is_valid_fasta()) {
        throw std::runtime_error {"Invalid FASTA " + fasta_path};
    }
    
//    fasta_index_path_ = fasta_path_;
//    fasta_index_path_.replace_extension("fai");
    
    fasta_index_path_ = fasta_path + ".fai";
    
    open_files();
}

Fasta::Fasta(std::string fasta_path, std::string fasta_index_path)
: fasta_path_ {std::move(fasta_path)}, fasta_index_path_ {std::move(fasta_index_path)}
{
    if (!is_valid_fasta()) {
        throw std::runtime_error {"Invalid FASTA " + fasta_path};
    }
    
    open_files();
}

bool Fasta::is_valid_fasta() const
{
//    auto extension = fasta_path_.extension().string();
//    if (!(extension == "fa" || extension == "fasta")) {
//        return false;
//    }
    return true; // TODO: could actually check valid fasta format
}

void Fasta::open_files()
{
    fasta_       = std::ifstream(fasta_path_);
    fasta_index_ = std::ifstream(fasta_index_path_);
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
    return static_cast<uint_fast32_t>(bioio::get_contig_size(fasta_index_, contig_name));
}

std::string Fasta::get_sequence(const GenomicRegion& a_region)
{
    if (is_in_cache(a_region)) {
        return region_cache_.at(a_region);
    } else {
        //TODO: this is a bit hacky/inefficient.. will change when bioio is better.
        auto contig_map = bioio::read_fasta_index(fasta_index_path_);
        auto index = contig_map[a_region.get_contig_name()];
        return bioio::read_fasta_contig(fasta_, index, a_region.get_begin(), size(a_region));
    }
}
