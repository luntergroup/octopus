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
#include <algorithm> // std::min, std::max
#include <stdexcept>

#include "comparable.h"

/*
    Represents a region of continuous sequence.
    begin and end positions are zero-indexed half open - [begin,end) - indexes.
 */
class SequenceRegion : public Comparable<SequenceRegion>
{
public:
    using SizeType       = std::uint_fast32_t;
    using DifferenceType = std::int_fast64_t;
    
    SequenceRegion() = default;
    explicit SequenceRegion(SizeType begin, SizeType end);
    
    SequenceRegion(const SequenceRegion&)            = default;
    SequenceRegion& operator=(const SequenceRegion&) = default;
    SequenceRegion(SequenceRegion&&)                 = default;
    SequenceRegion& operator=(SequenceRegion&&)      = default;
    
    SizeType get_begin() const noexcept;
    SizeType get_end() const noexcept;
    
private:
    SizeType begin_;
    SizeType end_;
};

inline SequenceRegion::SequenceRegion(SizeType begin, SizeType end)
:
begin_ {begin},
end_ {end}
{
    if (end < begin) throw std::runtime_error {"Invalid sequence region"};
}

inline SequenceRegion::SizeType SequenceRegion::get_begin() const noexcept
{
    return begin_;
}

inline SequenceRegion::SizeType SequenceRegion::get_end() const noexcept
{
    return end_;
}

inline SequenceRegion::SizeType size(const SequenceRegion& a_region)
{
    return a_region.get_end() - a_region.get_begin();
}

inline SequenceRegion::DifferenceType overlap_size(const SequenceRegion& lhs, const SequenceRegion& rhs)
{
    return static_cast<SequenceRegion::DifferenceType>(std::min(lhs.get_end(), rhs.get_end())) -
            static_cast<SequenceRegion::DifferenceType>(std::max(lhs.get_begin(), rhs.get_begin()));
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
    return (lhs.get_begin() == rhs.get_begin()) ? lhs.get_end() < rhs.get_end() : lhs.get_begin() < rhs.get_begin();
}

#endif /* defined(__Octopus__sequence_region__) */
