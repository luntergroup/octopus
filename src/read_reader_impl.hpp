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
#include <cstddef> // std::size_t
#include <cstdint> // std::uint_fast32_t
#include <unordered_map>

#include "genomic_region.hpp"

class AlignedRead;

class IReadReaderImpl
{
public:
    using SampleIdType       = std::string;
    using SizeType           = GenomicRegion::SizeType;
    using SampleIdToReadsMap = std::unordered_map<SampleIdType, std::vector<AlignedRead>>;
    
    virtual std::vector<SampleIdType> get_samples() = 0;
    virtual std::vector<std::string> get_read_groups_in_sample(const SampleIdType& sample) = 0;
    virtual std::size_t get_num_reads(const GenomicRegion& region) = 0;
    virtual SampleIdToReadsMap fetch_reads(const GenomicRegion& region) = 0;
    virtual unsigned get_num_reference_contigs() noexcept = 0;
    virtual std::vector<std::string> get_reference_contig_names() = 0;
    virtual SizeType get_reference_contig_size(const std::string& contig_name) = 0;
    virtual std::vector<GenomicRegion> get_possible_regions_in_file() = 0;
    
    virtual void open() = 0;
    virtual void close() = 0;
    virtual ~IReadReaderImpl() noexcept = default;
};

#endif
