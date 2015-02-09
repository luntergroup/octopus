//
//  sequence_region.h
//  Octopus
//
//  Created by Daniel Cooke on 09/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__sequence_region__
#define __Octopus__sequence_region__

#include "common.h"

/*
    Represents a region of continuous sequence.
    begin and end positions are zero-indexed half open - [begin,end) - indexes.
 */
struct SequenceRegion
{
    const int_fast32_t begin_pos, end_pos;
    
    SequenceRegion() = delete;
    SequenceRegion(int_fast32_t begin_pos, int_fast32_t end_pos);
};

int_fast32_t size(const SequenceRegion& a_region) noexcept;

inline bool operator==(const SequenceRegion& lhs, const SequenceRegion& rhs) noexcept;
inline bool operator!=(const SequenceRegion& lhs, const SequenceRegion& rhs) noexcept;
inline bool operator< (const SequenceRegion& lhs, const SequenceRegion& rhs) noexcept;
inline bool operator> (const SequenceRegion& lhs, const SequenceRegion& rhs) noexcept;
inline bool operator<=(const SequenceRegion& lhs, const SequenceRegion& rhs) noexcept;
inline bool operator>=(const SequenceRegion& lhs, const SequenceRegion& rhs) noexcept;

int_fast32_t overlap_size(const SequenceRegion& lhs, const SequenceRegion& rhs) noexcept;
bool overlaps(const SequenceRegion& lhs, const SequenceRegion& rhs) noexcept;

#endif /* defined(__Octopus__sequence_region__) */
