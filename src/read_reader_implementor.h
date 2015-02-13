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
#include <set>
#include <unordered_set>
#include <vector>
#include <cstdint>
#include <unordered_map>

#include "aligned_read.h"
#include "genomic_region.h"

using std::uint_fast32_t;

class IReadReaderImplementor
{
public:
    virtual std::set<AlignedRead> fetch_reads(const GenomicRegion& a_region) = 0;
    virtual uint_fast32_t get_num_reference_contigs() noexcept = 0;
    virtual std::unordered_set<std::string> get_reference_contig_names() = 0;
    virtual uint_fast32_t get_reference_contig_size(const std::string& contig_name) = 0;
    virtual ~IReadReaderImplementor() = default;
};

#endif
