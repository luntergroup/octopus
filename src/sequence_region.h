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
    if (end < begin) throw std::runtime_error {"Invalid sequence region: end < begin"};
}

inline SequenceRegion::SizeType SequenceRegion::get_begin() const noexcept
{
    return begin_;
}

inline SequenceRegion::SizeType SequenceRegion::get_end() const noexcept
{
    return end_;
}

inline bool empty(const SequenceRegion& a_region) noexcept
{
    return a_region.get_begin() == a_region.get_end();
}

inline SequenceRegion::SizeType size(const SequenceRegion& a_region) noexcept
{
    return a_region.get_end() - a_region.get_begin();
}

inline bool begins_equal(const SequenceRegion& lhs, const SequenceRegion& rhs) noexcept
{
    return lhs.get_begin() == rhs.get_begin();
}

inline bool ends_equal(const SequenceRegion& lhs, const SequenceRegion& rhs) noexcept
{
    return lhs.get_end() == rhs.get_end();
}

inline bool begins_before(const SequenceRegion& lhs, const SequenceRegion& rhs) noexcept
{
    return lhs.get_begin() < rhs.get_begin();
}

inline bool ends_before(const SequenceRegion& lhs, const SequenceRegion& rhs) noexcept
{
    return lhs.get_end() < rhs.get_end();
}

inline bool operator==(const SequenceRegion& lhs, const SequenceRegion& rhs) noexcept
{
    return begins_equal(lhs, rhs) && ends_equal(lhs, rhs);
}

inline bool operator<(const SequenceRegion& lhs, const SequenceRegion& rhs) noexcept
{
    return (begins_equal(lhs, rhs)) ? ends_before(lhs, rhs) : begins_before(lhs, rhs);
}

inline bool is_before(const SequenceRegion& lhs, const SequenceRegion& rhs) noexcept
{
    return (lhs == rhs) ? false : lhs.get_end() <= rhs.get_begin();
}

inline bool is_after(const SequenceRegion& lhs, const SequenceRegion& rhs) noexcept
{
    return  (lhs == rhs) ? false : rhs.get_end() <= lhs.get_begin();
}

inline bool are_adjacent(const SequenceRegion& lhs, const SequenceRegion& rhs) noexcept
{
    return lhs.get_begin() == rhs.get_end() || lhs.get_end() == rhs.get_begin();
}

inline SequenceRegion::DifferenceType overlap_size(const SequenceRegion& lhs, const SequenceRegion& rhs) noexcept
{
    return static_cast<SequenceRegion::DifferenceType>(std::min(lhs.get_end(), rhs.get_end())) -
            static_cast<SequenceRegion::DifferenceType>(std::max(lhs.get_begin(), rhs.get_begin()));
}

inline bool overlaps(const SequenceRegion& lhs, const SequenceRegion& rhs) noexcept
{
    auto num_bases_overlaped = overlap_size(lhs, rhs);
    return (num_bases_overlaped == 0) ? empty(lhs) || empty(rhs) : num_bases_overlaped > 0;
}

inline bool contains(const SequenceRegion& lhs, const SequenceRegion& rhs) noexcept
{
    return lhs.get_begin() <= rhs.get_begin() && rhs.get_end() <= lhs.get_end();
}

inline SequenceRegion get_overlapped(const SequenceRegion& lhs, const SequenceRegion& rhs) noexcept
{
    if (!overlaps(lhs, rhs)) {
        throw std::runtime_error {"Cannot get overlapped region between non overlapping regions"};
    }
    
    return SequenceRegion {std::max(lhs.get_begin(), rhs.get_begin()), std::min(lhs.get_end(), rhs.get_end())};
}

inline SequenceRegion get_intervening_region(const SequenceRegion& lhs, const SequenceRegion& rhs)
{
    if (begins_before(rhs, lhs) || overlaps(lhs, rhs)) {
        throw std::runtime_error {"Cannot get intervening region between overlapping regions"};
    }
    
    return SequenceRegion {lhs.get_end(), rhs.get_begin()};
}

inline SequenceRegion get_left_overhang(const SequenceRegion& lhs, const SequenceRegion& rhs) noexcept
{
    if (begins_before(rhs, lhs)) return SequenceRegion {lhs.get_begin(), lhs.get_begin()};
    
    return SequenceRegion {lhs.get_begin(), rhs.get_begin()};
}

inline SequenceRegion get_right_overhang(const SequenceRegion& lhs, const SequenceRegion& rhs) noexcept
{
    if (ends_before(lhs, rhs)) return SequenceRegion {lhs.get_end(), lhs.get_end()};
    
    return SequenceRegion {rhs.get_end(), lhs.get_end()};
}

#endif /* defined(__Octopus__sequence_region__) */
