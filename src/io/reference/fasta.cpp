// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "fasta.hpp"

#include <iostream>
#include <utility>
#include <exception>

#include <boost/filesystem/operations.hpp>

#include "basics/genomic_region.hpp"
#include "utils/sequence_utils.hpp"
#include "exceptions/missing_file_error.hpp"
#include "exceptions/missing_index_error.hpp"
#include "exceptions/malformed_file_error.hpp"
#include "exceptions/program_error.hpp"

namespace octopus { namespace io {

class MissingFasta : public MissingFileError
{
    std::string do_where() const override
    {
        return "Fasta";
    }
public:
    MissingFasta(Fasta::Path file) : MissingFileError {std::move(file), "fasta"} {}
};

class MalformedFasta : public MalformedFileError
{
    std::string do_where() const override
    {
        return "Fasta";
    }
public:
    MalformedFasta(Fasta::Path file) : MalformedFileError {std::move(file), "fasta"} {}
};

class MissingFastaIndex : public MissingIndexError
{
    std::string do_where() const override
    {
        return "Fasta";
    }
    
    std::string do_help() const override
    {
        return "ensure that a valid fasta index (.fai) exists in the same directory as the given "
        "fasta file. You can make one with the 'samtools faidx' command";
    }
public:
    MissingFastaIndex(Fasta::Path file) : MissingIndexError {std::move(file), "fasta"} {}
};

class MalformedFastaIndex : public MalformedFileError
{
    std::string do_where() const override
    {
        return "Fasta";
    }
public:
    MalformedFastaIndex(Fasta::Path file) : MalformedFileError {std::move(file), "fasta"} {}
};

Fasta::Fasta(Path fasta_path)
: Fasta {fasta_path, fasta_path.string() + ".fai", Options {}}
{}

Fasta::Fasta(Path fasta_path, Options options)
: Fasta {fasta_path, fasta_path.string() + ".fai", options}
{}

Fasta::Fasta(Path fasta_path, Path fasta_index_path)
: Fasta {std::move(fasta_path), std::move(fasta_index_path), Options {}}
{}

Fasta::Fasta(Path fasta_path, Path fasta_index_path, Options options)
: path_ {std::move(fasta_path)}
, index_path_ {std::move(fasta_index_path)}
, options_ {options}
{
    using boost::filesystem::exists;
    if (!exists(path_)) {
        throw MissingFasta {path_};
    }
    if (!is_valid_fasta()) {
        throw MalformedFasta {path_};
    }
    if (!exists(index_path_)) {
        index_path_ = path_;
        index_path_.replace_extension("fai");
        if (!exists(index_path_)) {
            throw MissingFastaIndex {path_};
        }
    }
    if (!is_valid_fasta_index()) {
        throw MalformedFastaIndex {index_path_};
    }
    fasta_       = std::ifstream(path_.string());
    fasta_index_ = bioio::read_fasta_index(index_path_.string());
}

Fasta::Fasta(const Fasta& other)
: path_ {other.path_}
, index_path_ {other.index_path_}
, fasta_ {path_.string()}
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

class BadReferenceRequestRegion : public ProgramError
{
    GenomicRegion region;
    
    std::string do_why() const override
    {
        return "Requested bad reference region " + to_string(region);
    }
    std::string do_help() const override
    {
        return "Send a debug report";
    }
    std::string do_where() const override
    {
        return "Fasta";
    }
public:
    BadReferenceRequestRegion(GenomicRegion region) : region {std::move(region)} {}
};

Fasta::GeneticSequence Fasta::do_fetch_sequence(const GenomicRegion& region) const
{
    try {
        auto result = bioio::read_fasta_contig(fasta_, fasta_index_.at(contig_name(region)),
                                               mapped_begin(region), size(region));
        if (is_capitalisation_requested()) {
            utils::capitalise(result);
        }
        if (result.size() < size(region)) {
            if (options_.base_fill_policy == Options::BaseFillPolicy::throw_exception) {
                throw BadReferenceRequestRegion {region};
            }
            if (options_.base_fill_policy == Options::BaseFillPolicy::fill_with_ns) {
                result.resize(size(region), 'N');
            }
        }
        return result;
    } catch (const std::ios::failure& e) {
        throw; // TODO: test to see if the file has disappeared, throw octopus error if so.
    } catch (const std::exception& e) {
        throw;
    } catch (...) {
        throw;
    }
}

bool Fasta::is_valid_fasta() const noexcept
{
    const auto extension = path_.extension().string();
    if (extension != ".fa" && extension != ".fasta") {
        return false;
    }
    return true; // TODO: could actually check valid fasta format
}

bool Fasta::is_valid_fasta_index() const noexcept
{
    const auto extension = index_path_.extension().string();
    if (extension != ".fai") {
        return false;
    }
    return true; // TODO: could actually check valid fasta format
}

bool Fasta::is_capitalisation_requested() const noexcept
{
    return options_.base_transform_policy == Options::BaseTransformPolicy::capitalise;
}

} // namespace io
} // namespace octopus
