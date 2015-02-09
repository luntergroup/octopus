//
//  genomic_region.h
//  Octopus
//
//  Created by Daniel Cooke on 09/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__genomic_region__
#define __Octopus__genomic_region__

#include <string>

#include "common.h"
#include "sequence_region.h"

/*
    Represents a continuous sequence region in a genome. The sequence
    name is the reference sequence name (usually a chromosome), and the
    begin and end positions are zero-indexed half open - [begin,end) - indexes.
 */
class GenomicRegion
{
public:
    GenomicRegion() = delete;
    GenomicRegion(std::string sequence_name, int_fast32_t begin, int_fast32_t end);
    //GenomicRegion(std::string the_region);
    
    const std::string& get_sequence_name() const noexcept;
    const SequenceRegion& get_sequence_region() const noexcept;
    int_fast32_t get_begin_pos() const noexcept;
    int_fast32_t get_end_pos() const noexcept;

private:
    std::string the_sequence_name_;
    SequenceRegion the_region_;
};

int_fast32_t size(const SequenceRegion& a_region) noexcept;

inline bool operator==(const GenomicRegion& lhs, const GenomicRegion& rhs) noexcept;
inline bool operator!=(const GenomicRegion& lhs, const GenomicRegion& rhs) noexcept;

std::string to_string(const GenomicRegion& a_region);

int_fast32_t overlap_size(const GenomicRegion& lhs, const GenomicRegion& rhs) noexcept;
bool overlaps(const GenomicRegion& lhs, const GenomicRegion& rhs) noexcept;

#endif /* defined(__Octopus__genomic_region__) */
