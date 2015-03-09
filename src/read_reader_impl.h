//
//  read_reader_implementor.h
//  Octopus
//
//  Created by Daniel Cooke on 12/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_read_reader_impl_h
#define Octopus_read_reader_impl_h

#include <string>
#include <vector>
#include <cstddef> // std::size_t
#include <cstdint> // std::uint_fast32_t
#include <unordered_map>

class AlignedRead;
class GenomicRegion;

class IReadReaderImpl
{
public:
    using SampleIdType       = std::string;
    using SizeType           = std::uint_fast32_t;
    using SampleIdToReadsMap = std::unordered_map<SampleIdType, std::vector<AlignedRead>>;
    
    virtual std::vector<SampleIdType> get_sample_ids() = 0;
    virtual std::vector<std::string> get_read_groups_in_sample(const SampleIdType& a_sample_id) = 0;
    virtual std::size_t get_num_reads(const GenomicRegion& a_region) = 0;
    virtual SampleIdToReadsMap fetch_reads(const GenomicRegion& a_region) = 0;
    virtual SizeType get_num_reference_contigs() noexcept = 0;
    virtual std::vector<std::string> get_reference_contig_names() = 0;
    virtual SizeType get_reference_contig_size(const std::string& contig_name) = 0;
    virtual std::vector<GenomicRegion> get_regions_in_file() = 0;
    virtual void close() = 0; // may throw
    virtual ~IReadReaderImpl() noexcept = default;
};

#endif
