//
//  read_reader_implementor.hpp
//  Octopus
//
//  Created by Daniel Cooke on 12/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_read_reader_impl_hpp
#define Octopus_read_reader_impl_hpp

#include <string>
#include <vector>
#include <cstddef>
#include <unordered_map>

#include "genomic_region.hpp"

#include "common.hpp"

class AlignedRead;

class IReadReaderImpl
{
public:
    using SampleIdType  = std::string;
    using SizeType      = GenomicRegion::SizeType;
    using Reads         = Octopus::ReadContainer;
    using SampleReadMap = Octopus::ReadMap;
    
    virtual ~IReadReaderImpl() noexcept = default;
    
    virtual void open() = 0;
    virtual bool is_open() const noexcept= 0;
    virtual void close() = 0;
    
    virtual std::vector<SampleIdType> get_samples() = 0;
    virtual std::vector<std::string> get_read_groups_in_sample(const SampleIdType& sample) = 0;
    
    virtual size_t count_reads(const GenomicRegion& region) = 0;
    virtual size_t count_reads(const SampleIdType& sample, const GenomicRegion& region) = 0;
    virtual GenomicRegion find_covered_subregion(const GenomicRegion& region, size_t target_coverage) = 0;
    virtual SampleReadMap fetch_reads(const GenomicRegion& region) = 0;
    virtual Reads fetch_reads(const SampleIdType& sample, const GenomicRegion& region) = 0;
    
    virtual unsigned get_num_reference_contigs() noexcept = 0;
    virtual std::vector<std::string> get_reference_contig_names() = 0;
    virtual SizeType get_reference_contig_size(const std::string& contig_name) = 0;
    virtual std::vector<GenomicRegion> get_possible_regions_in_file() = 0;
};

#endif
