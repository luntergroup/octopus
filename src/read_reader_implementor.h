//
//  read_reader_implementor.h
//  Octopus
//
//  Created by Daniel Cooke on 12/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_read_reader_implementor_h
#define Octopus_read_reader_implementor_h

#include <string>
#include <vector>
#include <cstdint>
#include <unordered_map>

#include "aligned_read.h"
#include "genomic_region.h"

using std::uint_fast32_t;

class IReadReaderImplementor
{
public:
    using SampleIdToReadsMap = std::unordered_map<std::string, std::vector<AlignedRead>>;
    
    virtual std::vector<std::string> get_sample_ids() = 0;
    virtual std::vector<std::string> get_read_groups_in_sample(const std::string& a_sample_id) = 0;
    virtual SampleIdToReadsMap fetch_reads(const GenomicRegion& a_region) = 0;
    virtual uint_fast32_t get_num_reference_contigs() noexcept = 0;
    virtual std::vector<std::string> get_reference_contig_names() = 0;
    virtual uint_fast32_t get_reference_contig_size(const std::string& contig_name) = 0;
    virtual std::vector<GenomicRegion> get_regions_in_file() = 0;
    virtual void close() = 0; // may throw
    virtual ~IReadReaderImplementor() noexcept = default;
};

#endif
