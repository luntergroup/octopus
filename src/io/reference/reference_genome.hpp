// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef reference_genome_hpp
#define reference_genome_hpp

#include <string>
#include <vector>
#include <unordered_map>
#include <cstddef>
#include <memory>

#include <boost/filesystem/path.hpp>

#include "basics/genomic_region.hpp"
#include "utils/memory_footprint.hpp"
#include "reference_reader.hpp"

namespace octopus {

class ReferenceGenome
{
public:
    using ContigName      = io::ReferenceReader::ContigName;
    using GeneticSequence = io::ReferenceReader::GeneticSequence;
    
    ReferenceGenome() = delete;
    
    ReferenceGenome(std::unique_ptr<io::ReferenceReader> impl);
    
    ReferenceGenome(const ReferenceGenome&);
    ReferenceGenome& operator=(ReferenceGenome);
    ReferenceGenome(ReferenceGenome&&)            = default;
    ReferenceGenome& operator=(ReferenceGenome&&) = default;
    
    ~ReferenceGenome() = default;
    
    const std::string& name() const;
    
    bool has_contig(const ContigName& contig) const noexcept;
    std::size_t num_contigs() const noexcept;
    std::vector<ContigName> contig_names() const;
    ContigRegion::Size contig_size(const ContigName& contig) const;
    GenomicRegion contig_region(const ContigName& contig) const;
    
    bool contains(const GenomicRegion& region) const noexcept;
    
    GeneticSequence fetch_sequence(const GenomicRegion& region) const;
    
private:
    std::unique_ptr<io::ReferenceReader> impl_;
    std::string name_;
    std::unordered_map<ContigName, ContigRegion::Size> contig_sizes_;
    std::vector<ContigName> ordered_contigs_;
};

// non-member functions

ReferenceGenome make_reference(boost::filesystem::path reference_path,
                               MemoryFootprint max_cache_size = 0,
                               bool is_threaded = false,
                               bool capitalise_bases = true);

std::vector<GenomicRegion> get_all_contig_regions(const ReferenceGenome& reference);

GenomicRegion::Position calculate_genome_size(const ReferenceGenome& reference);

} // namespace octopus

#endif
