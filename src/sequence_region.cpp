//
//  sequence_region.cpp
//  Octopus
//
//  Created by Daniel Cooke on 09/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "sequence_region.h"

#include <algorithm>

SequenceRegion::SequenceRegion(int_fast32_t begin_pos, int_fast32_t end_pos)
: begin_pos {begin_pos}, end_pos {end_pos}
{}

int_fast32_t size(const SequenceRegion& a_region) noexcept
{
    return a_region.end_pos - a_region.begin_pos;
}

inline bool operator==(const SequenceRegion& lhs, const SequenceRegion& rhs) noexcept
{
    return lhs.begin_pos == rhs.begin_pos && lhs.end_pos == rhs.end_pos;
}
inline bool operator!=(const SequenceRegion& lhs, const SequenceRegion& rhs) noexcept {return !operator==(lhs, rhs);}
inline bool operator< (const SequenceRegion& lhs, const SequenceRegion& rhs) noexcept
{
    return lhs.begin_pos < rhs.begin_pos;
}
inline bool operator> (const SequenceRegion& lhs, const SequenceRegion& rhs) noexcept {return operator<(rhs,lhs);}
inline bool operator<=(const SequenceRegion& lhs, const SequenceRegion& rhs) noexcept {return !operator>(lhs,rhs);}
inline bool operator>=(const SequenceRegion& lhs, const SequenceRegion& rhs) noexcept {return !operator<(lhs,rhs);}

int_fast32_t overlap_size(const SequenceRegion& lhs, const SequenceRegion& rhs) noexcept
{
    return std::min(lhs.end_pos, rhs.end_pos) - std::max(lhs.begin_pos, rhs.begin_pos);
}

bool overlaps(const SequenceRegion& lhs, const SequenceRegion& rhs) noexcept
{
    return overlap_size(lhs, rhs) > 0;
}
