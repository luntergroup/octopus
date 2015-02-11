//
//  aligned_read.h
//  Octopus
//
//  Created by Daniel Cooke on 11/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__aligned_read__
#define __Octopus__aligned_read__

#include <string>
#include <cstdint>

#include "genomic_region.h"

using std::uint_fast32_t;

class AlignedRead
{
public:
    std::string get_sequence();
    std::string get_qualities();
    uint_fast32_t get_sequence_size();
    std::string get_cigar_string();
    GenomicRegion get_region();
    std::string get_contig_name();
    uint_fast32_t get_begin();
    uint_fast32_t get_end();
    uint_fast32_t get_insert_size();
    std::string get_mate_contig_name();
    uint_fast32_t get_mate_begin();
    
private:
    
};

#endif /* defined(__Octopus__aligned_read__) */
