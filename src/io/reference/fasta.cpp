//
//  fasta.cpp
//  Octopus
//
//  Created by Daniel Cooke on 10/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "fasta.hpp"

#include <iostream>
#include <sstream>
#include <utility>

#include <boost/filesystem/operations.hpp>

#include <basics/genomic_region.hpp>

namespace octopus { namespace io {

Fasta::Fasta(Path fasta_path)
:
Fasta {fasta_path, fasta_path.string() + ".fai"}
{}

Fasta::Fasta(Path fasta_path, Path fasta_index_path)
:
path_ {std::move(fasta_path)},
index_path_ {std::move(fasta_index_path)}
{
    using boost::filesystem::exists;
    
    if (!exists(path_)) {
        std::ostringstream ss {};
        ss << "Fasta: given fasta path ";
        ss << fasta_path;
        ss << " does not exist";
        throw std::runtime_error {ss.str()};
    }
    
    if (!exists(index_path_)) {
        index_path_ = path_;
        index_path_.replace_extension("fai");
        
        if (!exists(index_path_)) {
            throw MissingFastaIndex {path_, index_path_};
        }
    }
    
    if (is_valid()) {
        fasta_       = std::ifstream(path_.string());
        fasta_index_ = bioio::read_fasta_index(index_path_.string());
    }
}

Fasta::Fasta(const Fasta& other)
:
path_ {other.path_},
index_path_ {other.index_path_},
fasta_ {path_.string()}
{}

Fasta& Fasta::operator=(Fasta other)
{
    using std::swap;
    swap(path_, other.path_);
    swap(index_path_, other.index_path_);
    swap(fasta_, other.fasta_);
    return *this;
}

// virtual private methods

std::unique_ptr<ReferenceReader> Fasta::do_clone() const
{
    return std::make_unique<Fasta>(*this);
}

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
    return path_.stem().string();
}

std::vector<Fasta::ContigName> Fasta::do_fetch_contig_names() const
{
    return bioio::read_fasta_index_contig_names(index_path_.string());
}

Fasta::GenomicSize Fasta::do_fetch_contig_size(const ContigName& contig) const
{
    if (fasta_index_.count(contig) == 0) {
        throw std::runtime_error {"contig \"" + contig +
            "\" not found in fasta index \"" + index_path_.string() + "\""};
    }
    return static_cast<GenomicSize>(fasta_index_.at(contig).length);
}

Fasta::GeneticSequence Fasta::do_fetch_sequence(const GenomicRegion& region) const
{
    return bioio::read_fasta_contig(fasta_, fasta_index_.at(contig_name(region)),
                                    mapped_begin(region), size(region));
}

bool Fasta::is_valid() const noexcept
{
    const auto extension = path_.extension().string();
    
    if (extension != ".fa" && extension != ".fasta") {
        return false;
    }
    
    return true; // TODO: could actually check valid fasta format
}

// MissingFastaIndex

Fasta::MissingFastaIndex::MissingFastaIndex(Path fasta_path, Path expected_index_path)
:
std::runtime_error {"MissingFastaIndex"},
fasta_path_ {std::move(fasta_path)},
expected_index_path_ {expected_index_path},
msg_ {}
{}

const Fasta::MissingFastaIndex::Path& Fasta::MissingFastaIndex::fasta_path() const noexcept
{
    return fasta_path_;
}

const Fasta::MissingFastaIndex::Path& Fasta::MissingFastaIndex::expected_index_path() const noexcept
{
    return expected_index_path_;
}

const char* Fasta::MissingFastaIndex::what() const noexcept
{
    std::ostringstream ss {};
    ss << runtime_error::what() << ": expected index " << expected_index_path_ << " for FASTA file "
        << fasta_path_ << " but it is missing";
    msg_ = ss.str();
    return msg_.c_str();
}

} // namespace io
} // namespace octopus
