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
    using SizeType     = GenomicRegion::SizeType;
    using SequenceType = std::string;
    
    std::string get_reference_name() const;
    std::vector<std::string> get_contig_names();
    SizeType get_contig_size(const std::string& contig_name);
    SequenceType get_sequence(const GenomicRegion& region);
    
    virtual ~ReferenceGenomeImpl() noexcept = default;
    
private:
    virtual std::string do_get_reference_name() const = 0;
    virtual std::vector<std::string> do_get_contig_names() = 0;
    virtual SizeType do_get_contig_size(const std::string& contig_name) = 0;
    virtual SequenceType do_get_sequence(const GenomicRegion& region) = 0;
};

inline std::string ReferenceGenomeImpl::get_reference_name() const
{
    return do_get_reference_name();
}

inline std::vector<std::string> ReferenceGenomeImpl::get_contig_names()
{
    return do_get_contig_names();
}

inline ReferenceGenomeImpl::SizeType ReferenceGenomeImpl::get_contig_size(const std::string& contig_name)
{
    return do_get_contig_size(contig_name);
}

inline ReferenceGenomeImpl::SequenceType ReferenceGenomeImpl::get_sequence(const GenomicRegion& region)
{
    return do_get_sequence(region);
}

#endif
