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

using std::uint_fast32_t;
using std::int_fast64_t;

/*
    Represents a region of continuous sequence.
    begin and end positions are zero-indexed half open - [begin,end) - indexes.
 */
class SequenceRegion
{
public:
    SequenceRegion() = delete;
    SequenceRegion(uint_fast32_t begin_pos, uint_fast32_t end_pos);
    
    uint_fast32_t get_begin_pos() const noexcept;
    uint_fast32_t get_end_pos() const noexcept;
    
private:
    const uint_fast32_t begin_pos_;
    const uint_fast32_t end_pos_;
};

inline SequenceRegion::SequenceRegion(uint_fast32_t begin_pos, uint_fast32_t end_pos)
: begin_pos_ {begin_pos}, end_pos_ {end_pos}
{}

inline uint_fast32_t SequenceRegion::get_begin_pos() const noexcept
{
    return begin_pos_;
}

inline uint_fast32_t SequenceRegion::get_end_pos() const noexcept
{
    return end_pos_;
}

inline uint_fast32_t size(const SequenceRegion& a_region)
{
    return a_region.get_end_pos() - a_region.get_begin_pos();
}

inline int_fast64_t overlap_size(const SequenceRegion& lhs, const SequenceRegion& rhs)
{
    return static_cast<int_fast64_t>(std::min(lhs.get_end_pos(), rhs.get_end_pos())) -
            static_cast<int_fast64_t>(std::max(lhs.get_begin_pos(), rhs.get_begin_pos()));
}

inline bool overlaps(const SequenceRegion& lhs, const SequenceRegion& rhs)
{
    return overlap_size(lhs, rhs) > 0;
}

inline bool operator==(const SequenceRegion& lhs, const SequenceRegion& rhs)
{
    return lhs.get_begin_pos() == rhs.get_begin_pos() && lhs.get_end_pos() == rhs.get_end_pos();
}
inline bool operator< (const SequenceRegion& lhs, const SequenceRegion& rhs)
{
    return lhs.get_begin_pos() < rhs.get_begin_pos();
}
inline bool operator!=(const SequenceRegion& lhs, const SequenceRegion& rhs){return !operator==(lhs, rhs);}
inline bool operator> (const SequenceRegion& lhs, const SequenceRegion& rhs){return operator<(rhs,lhs);}
inline bool operator<=(const SequenceRegion& lhs, const SequenceRegion& rhs){return !operator>(lhs,rhs);}
inline bool operator>=(const SequenceRegion& lhs, const SequenceRegion& rhs){return !operator<(lhs,rhs);}

#endif /* defined(__Octopus__sequence_region__) */
