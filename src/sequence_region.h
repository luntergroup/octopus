//
//  sequence_region.h
//  Octopus
//
//  Created by Daniel Cooke on 09/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__sequence_region__
#define __Octopus__sequence_region__

#include <cstdint>

#include "comparable.h"

using std::uint_fast32_t;
using std::int_fast64_t;

/*
    Represents a region of continuous sequence.
    begin and end positions are zero-indexed half open - [begin,end) - indexes.
 */
class SequenceRegion : Comparable<SequenceRegion>
{
public:
    SequenceRegion() = delete;
    SequenceRegion(uint_fast32_t begin, uint_fast32_t end);
    
    SequenceRegion(const SequenceRegion&)            = default;
    SequenceRegion& operator=(const SequenceRegion&) = default;
    SequenceRegion(SequenceRegion&&)                 = default;
    SequenceRegion& operator=(SequenceRegion&&)      = default;
    
    uint_fast32_t get_begin() const noexcept;
    uint_fast32_t get_end() const noexcept;
    
private:
    uint_fast32_t begin_;
    uint_fast32_t end_;
};

inline SequenceRegion::SequenceRegion(uint_fast32_t begin, uint_fast32_t end)
:begin_ {begin},
 end_ {end}
{}

inline uint_fast32_t SequenceRegion::get_begin() const noexcept
{
    return begin_;
}

inline uint_fast32_t SequenceRegion::get_end() const noexcept
{
    return end_;
}

inline uint_fast32_t size(const SequenceRegion& a_region)
{
    return a_region.get_end() - a_region.get_begin();
}

inline int_fast64_t overlap_size(const SequenceRegion& lhs, const SequenceRegion& rhs)
{
    return static_cast<int_fast64_t>(std::min(lhs.get_end(), rhs.get_end())) -
            static_cast<int_fast64_t>(std::max(lhs.get_begin(), rhs.get_begin()));
}

inline bool overlaps(const SequenceRegion& lhs, const SequenceRegion& rhs)
{
    return overlap_size(lhs, rhs) > 0;
}

inline bool operator==(const SequenceRegion& lhs, const SequenceRegion& rhs)
{
    return lhs.get_begin() == rhs.get_begin() && lhs.get_end() == rhs.get_end();
}

inline bool operator<(const SequenceRegion& lhs, const SequenceRegion& rhs)
{
    return lhs.get_begin() < rhs.get_begin();
}

#endif /* defined(__Octopus__sequence_region__) */
