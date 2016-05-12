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
#include <utility>

#include "genomic_region.hpp"
#include "aligned_read.hpp"

class IReadReaderImpl
{
public:
    using SampleIdType  = std::string;
    using SizeType      = GenomicRegion::SizeType;
    using ReadContainer = std::vector<AlignedRead>;
    using SampleReadMap = std::unordered_map<SampleIdType, ReadContainer>;
    
    using CoveragePair = std::pair<GenomicRegion, std::vector<unsigned>>;
    
    virtual ~IReadReaderImpl() noexcept = default;
    
    virtual bool is_open() const noexcept = 0;
    virtual void open() = 0;
    virtual void close() = 0;
    
    virtual std::vector<SampleIdType> extract_samples() = 0;
    virtual std::vector<std::string> extract_read_groups_in_sample(const SampleIdType& sample) = 0;
    
    virtual bool has_contig_reads(const GenomicRegion::ContigNameType& contig) = 0;
    
    virtual std::size_t count_reads(const GenomicRegion& region) = 0;
    virtual std::size_t count_reads(const SampleIdType& sample, const GenomicRegion& region) = 0;
    
    virtual CoveragePair find_covered_subregion(const GenomicRegion& region, std::size_t max_coverage) = 0;
    virtual CoveragePair find_covered_subregion(const SampleIdType& sample, const GenomicRegion& region,
                                                std::size_t max_coverage) = 0;
    virtual CoveragePair find_covered_subregion(const std::vector<SampleIdType>& samples, const GenomicRegion& region,
                                                std::size_t max_coverage) = 0;
    
    virtual SampleReadMap fetch_reads(const GenomicRegion& region) = 0;
    virtual ReadContainer fetch_reads(const SampleIdType& sample, const GenomicRegion& region) = 0;
    virtual SampleReadMap fetch_reads(const std::vector<SampleIdType>& samples, const GenomicRegion& region) = 0;
    
    virtual unsigned count_reference_contigs() = 0;
    virtual std::vector<std::string> extract_reference_contig_names() = 0;
    virtual SizeType get_reference_contig_size(const std::string& contig_name) = 0;
    virtual std::vector<GenomicRegion> extract_possible_regions_in_file() = 0;
};

#endif
