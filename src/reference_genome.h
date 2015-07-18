//
//  reference_genome.h
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
#include <cstddef> // std::size_t
#include <memory>  // std::unique_ptr

#include "genomic_region.h"
#include "reference_genome_impl.h"

//TODO: make sure this is thread-safe

class ReferenceGenome
{
public:
    using SizeType     = IReferenceGenomeImpl::SizeType;
    using SequenceType = IReferenceGenomeImpl::SequenceType;
    
    ReferenceGenome() = delete;
    explicit ReferenceGenome(std::unique_ptr<IReferenceGenomeImpl> the_reference_impl);
    
    ReferenceGenome(const ReferenceGenome&)            = delete;
    ReferenceGenome& operator=(const ReferenceGenome&) = delete;
    ReferenceGenome(ReferenceGenome&&)                 = default;
    ReferenceGenome& operator=(ReferenceGenome&&)      = default;
    
    const std::string& get_name() const;
    bool has_contig(const std::string& contig_name) const noexcept;
    std::size_t num_contigs() const noexcept;
    const std::vector<std::string>& get_contig_names() const noexcept;
    SizeType get_contig_size(const std::string& contig_name) const;
    SizeType get_contig_size(const GenomicRegion& a_region) const;
    GenomicRegion get_contig_region(const std::string& contig_name) const;
    bool contains_region(const GenomicRegion& a_region) const noexcept;
    SequenceType get_sequence(const GenomicRegion& a_region);
    
private:
    std::unique_ptr<IReferenceGenomeImpl> the_reference_impl_;
    std::string name_;
    std::vector<std::string> contig_names_;
    std::unordered_map<std::string, SizeType> contig_sizes_;
};

// non-member functions

std::vector<GenomicRegion> get_all_contig_regions(const ReferenceGenome& the_reference);

// Requires reference access to get contig sizes for partially specified regions (e.g. "4")
GenomicRegion parse_region(const std::string& a_region, const ReferenceGenome& the_reference);

#endif /* defined(__Octopus__reference_genome__) */
