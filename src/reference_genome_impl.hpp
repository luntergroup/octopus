//
//  reference_genome_implementor.hpp
//  Octopus
//
//  Created by Daniel Cooke on 10/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_reference_genome_impl_hpp
#define Octopus_reference_genome_impl_hpp

#include <string>
#include <vector>
#include <cstdint>

#include "genomic_region.hpp"

class ReferenceGenomeImpl
{
public:
    using ContigName      = GenomicRegion::ContigName;
    using GenomicSize     = GenomicRegion::Size;
    using GeneticSequence = std::string;
    
    bool is_open() const noexcept;
    
    std::string fetch_reference_name() const;
    std::vector<ContigName> fetch_contig_names() const;
    GenomicSize fetch_contig_size(const ContigName& contig) const;
    GeneticSequence fetch_sequence(const GenomicRegion& region) const;
    
    virtual ~ReferenceGenomeImpl() noexcept = default;
    
private:
    virtual bool do_is_open() const noexcept = 0;
    virtual std::string do_fetch_reference_name() const = 0;
    virtual std::vector<ContigName> do_fetch_contig_names() const = 0;
    virtual GenomicSize do_fetch_contig_size(const ContigName& contig) const = 0;
    virtual GeneticSequence do_fetch_sequence(const GenomicRegion& region) const = 0;
};

inline bool ReferenceGenomeImpl::is_open() const noexcept
{
    return do_is_open();
}

inline std::string ReferenceGenomeImpl::fetch_reference_name() const
{
    return do_fetch_reference_name();
}

inline std::vector<ReferenceGenomeImpl::ContigName> ReferenceGenomeImpl::fetch_contig_names() const
{
    return do_fetch_contig_names();
}

inline ReferenceGenomeImpl::GenomicSize ReferenceGenomeImpl::fetch_contig_size(const ContigName& contig) const
{
    return do_fetch_contig_size(contig);
}

inline ReferenceGenomeImpl::GeneticSequence ReferenceGenomeImpl::fetch_sequence(const GenomicRegion& region) const
{
    return do_fetch_sequence(region);
}

#endif
