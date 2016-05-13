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
    using ContigNameType = GenomicRegion::ContigNameType;
    using SizeType       = GenomicRegion::SizeType;
    using SequenceType   = std::string;
    
    bool is_open() const noexcept;
    
    std::string fetch_reference_name() const;
    std::vector<ContigNameType> fetch_contig_names() const;
    SizeType fetch_contig_size(const ContigNameType& contig) const;
    SequenceType fetch_sequence(const GenomicRegion& region) const;
    
    virtual ~ReferenceGenomeImpl() noexcept = default;
    
private:
    virtual bool do_is_open() const noexcept = 0;
    virtual std::string do_fetch_reference_name() const = 0;
    virtual std::vector<ContigNameType> do_fetch_contig_names() const = 0;
    virtual SizeType do_fetch_contig_size(const ContigNameType& contig) const = 0;
    virtual SequenceType do_fetch_sequence(const GenomicRegion& region) const = 0;
};

inline bool ReferenceGenomeImpl::is_open() const noexcept
{
    return do_is_open();
}

inline std::string ReferenceGenomeImpl::fetch_reference_name() const
{
    return do_fetch_reference_name();
}

inline std::vector<ReferenceGenomeImpl::ContigNameType> ReferenceGenomeImpl::fetch_contig_names() const
{
    return do_fetch_contig_names();
}

inline ReferenceGenomeImpl::SizeType ReferenceGenomeImpl::fetch_contig_size(const ContigNameType& contig) const
{
    return do_fetch_contig_size(contig);
}

inline ReferenceGenomeImpl::SequenceType ReferenceGenomeImpl::fetch_sequence(const GenomicRegion& region) const
{
    return do_fetch_sequence(region);
}

#endif
