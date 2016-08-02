//
//  reference_genome.hpp
//  Octopus
//
//  Created by Daniel Cooke on 10/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__reference_genome__
#define __Octopus__reference_genome__

#include <string>
#include <vector>
#include <unordered_map>
#include <cstddef>
#include <memory>

#include <boost/filesystem/path.hpp>

#include <basics/genomic_region.hpp>

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
    ReferenceGenome(ReferenceGenome&&)                 = default;
    ReferenceGenome& operator=(ReferenceGenome&&)      = default;
    
    ~ReferenceGenome() = default;
    
    bool is_good() const noexcept;
    
    const std::string& name() const;
    
    bool has_contig(const ContigName& contig) const noexcept;
    
    bool contains_region(const GenomicRegion& region) const noexcept;
    
    std::size_t num_contigs() const noexcept;
    
    const std::vector<ContigName>& contig_names() const noexcept;
    
    ContigRegion::Size contig_size(const ContigName& contig) const;
    
    ContigRegion::Size contig_size(const GenomicRegion& region) const;
    
    GenomicRegion contig_region(const ContigName& contig) const;
    
    GeneticSequence fetch_sequence(const GenomicRegion& region) const;
    
private:
    std::unique_ptr<io::ReferenceReader> impl_;
    
    std::string name_;
    std::vector<ContigName> contig_names_;
    std::unordered_map<ContigName, ContigRegion::Size> contig_sizes_;
};

// non-member functions

ReferenceGenome make_reference(boost::filesystem::path reference_path,
                               std::size_t max_cached_bases = 0,
                               bool is_threaded = false);

std::vector<GenomicRegion> get_all_contig_regions(const ReferenceGenome& reference);

GenomicRegion::Position calculate_genome_size(const ReferenceGenome& reference);

// Requires reference access to get contig sizes for partially specified regions (e.g. "4")
GenomicRegion parse_region(std::string region, const ReferenceGenome& reference);

} // namespace octopus

#endif /* defined(__Octopus__reference_genome__) */
