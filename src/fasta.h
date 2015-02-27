//
//  fasta.h
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
#include <exception>
#include <fstream>
#include <unordered_map>
#include <boost/filesystem/path.hpp>
#include <boost/filesystem/operations.hpp>

#include "reference_genome_implementor.h"
#include "genomic_region.h"
#include "bioio.h"

namespace fs = boost::filesystem;

class Fasta : public IReferenceGenomeImplementor
{
public:
    Fasta() = delete;
    Fasta(std::string fasta_path);
    Fasta(std::string fasta_path, std::string fasta_index_path);
    ~Fasta() noexcept override = default;
    
    Fasta(const Fasta&)            = default;
    Fasta& operator=(const Fasta&) = default;
    Fasta(Fasta&&)                 = default;
    Fasta& operator=(Fasta&&)      = default;
    
    std::string get_reference_name() override;
    std::vector<std::string> get_contig_names() override;
    std::uint_fast32_t get_contig_size(std::string contig_name) override;
    std::string get_sequence(const GenomicRegion& a_region) override;

private:
    fs::path fasta_path_;
    fs::path fasta_index_path_;
    std::ifstream fasta_;
    std::unordered_map<std::string, bioio::FastaIndex> fasta_contig_indices_;
    std::unordered_map<GenomicRegion, std::string> region_cache_; // TODO: is this useful?
    
    bool is_valid_fasta() const;
    bool is_in_cache(const GenomicRegion& a_region) const noexcept;
};

inline bool Fasta::is_in_cache(const GenomicRegion& a_region) const noexcept
{
    return region_cache_.count(a_region) > 0;
}

#endif /* defined(__Octopus__fasta__) */
