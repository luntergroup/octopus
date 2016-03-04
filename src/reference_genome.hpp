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

#include "genomic_region.hpp"
#include "reference_genome_impl.hpp"

class ReferenceGenome
{
public:
    using ContigNameType = ReferenceGenomeImpl::ContigNameType;
    using SizeType       = ReferenceGenomeImpl::SizeType;
    using SequenceType   = ReferenceGenomeImpl::SequenceType;
    
    ReferenceGenome() = delete;
    explicit ReferenceGenome(std::unique_ptr<ReferenceGenomeImpl> impl);
    ~ReferenceGenome() = default;
    
    ReferenceGenome(const ReferenceGenome&)            = delete;
    ReferenceGenome& operator=(const ReferenceGenome&) = delete;
    ReferenceGenome(ReferenceGenome&&)                 = default;
    ReferenceGenome& operator=(ReferenceGenome&&)      = default;
    
    bool is_good() const noexcept;
    
    const std::string& get_name() const;
    
    bool has_contig(const ContigNameType& contig) const noexcept;
    bool contains_region(const GenomicRegion& region) const noexcept;
    std::size_t num_contigs() const noexcept;
    const std::vector<ContigNameType>& get_contig_names() const noexcept;
    SizeType get_contig_size(const ContigNameType& contig) const;
    SizeType get_contig_size(const GenomicRegion& region) const;
    GenomicRegion get_contig_region(const ContigNameType& contig) const;
    
    SequenceType get_sequence(const GenomicRegion& region) const;
    
private:
    std::unique_ptr<ReferenceGenomeImpl> impl_;
    
    std::string name_;
    std::vector<ContigNameType> contig_names_;
    std::unordered_map<ContigNameType, SizeType> contig_sizes_;
};

// non-member functions

ReferenceGenome make_reference(boost::filesystem::path reference_path,
                               std::size_t max_base_pair_cache = 0,
                               bool is_threaded = false);

std::vector<GenomicRegion> get_all_contig_regions(const ReferenceGenome& reference);

GenomicRegion::SizeType calculate_genome_size(const ReferenceGenome& reference);

// Requires reference access to get contig sizes for partially specified regions (e.g. "4")
GenomicRegion parse_region(std::string region, const ReferenceGenome& reference);

#endif /* defined(__Octopus__reference_genome__) */
