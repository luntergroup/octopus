//
//  contig_region.hpp
//  Octopus
//
//  Created by Daniel Cooke on 09/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__contig_region__
#define __Octopus__contig_region__

#include <ostream>
#include <cstdint>
#include <algorithm>
#include <stdexcept>

#include <boost/functional/hash.hpp>

#include "comparable.hpp"

/*
    Represents a region of continuous sequence.
    begin and end positions are zero-indexed half open indices [begin,end).
*/
class ContigRegion : public Comparable<ContigRegion>
{
public:
    using SizeType       = std::uint_fast32_t;
    using DifferenceType = std::int_fast64_t;
    
    ContigRegion() = default;
    explicit ContigRegion(SizeType begin, SizeType end);
    ~ContigRegion() = default;
    
    ContigRegion(const ContigRegion&)            = default;
    ContigRegion& operator=(const ContigRegion&) = default;
    ContigRegion(ContigRegion&&)                 = default;
    ContigRegion& operator=(ContigRegion&&)      = default;
    
    SizeType get_begin() const noexcept;
    SizeType get_end() const noexcept;
    
private:
    SizeType begin_, end_;
};

// public member methods

inline ContigRegion::ContigRegion(SizeType begin, SizeType end)
:
begin_ {begin},
end_ {end}
{
    if (end < begin) throw std::runtime_error {"ContigRegion constructed with end < begin"};
}

inline ContigRegion::SizeType ContigRegion::get_begin() const noexcept
{
    return begin_;
}

inline ContigRegion::SizeType ContigRegion::get_end() const noexcept
{
    return end_;
}

// non-member methods

inline ContigRegion::SizeType get_begin(const ContigRegion& region) noexcept
{
    return region.get_begin();
}

inline ContigRegion::SizeType get_end(const ContigRegion& region) noexcept
{
    return region.get_end();
}

inline bool empty(const ContigRegion& region) noexcept
{
    return region.get_begin() == region.get_end();
}

inline ContigRegion::SizeType size(const ContigRegion& region) noexcept
{
    return region.get_end() - region.get_begin();
}

inline bool begins_equal(const ContigRegion& lhs, const ContigRegion& rhs) noexcept
{
    return lhs.get_begin() == rhs.get_begin();
}

inline bool ends_equal(const ContigRegion& lhs, const ContigRegion& rhs) noexcept
{
    return lhs.get_end() == rhs.get_end();
}

inline bool begins_before(const ContigRegion& lhs, const ContigRegion& rhs) noexcept
{
    return lhs.get_begin() < rhs.get_begin();
}

inline bool ends_before(const ContigRegion& lhs, const ContigRegion& rhs) noexcept
{
    return lhs.get_end() < rhs.get_end();
}

inline bool operator==(const ContigRegion& lhs, const ContigRegion& rhs) noexcept
{
    return begins_equal(lhs, rhs) && ends_equal(lhs, rhs);
}

inline bool operator<(const ContigRegion& lhs, const ContigRegion& rhs) noexcept
{
    return begins_before(lhs, rhs) || (begins_equal(lhs, rhs) && ends_before(lhs, rhs));
}

inline bool is_before(const ContigRegion& lhs, const ContigRegion& rhs) noexcept
{
    return lhs.get_end() <= rhs.get_begin() && lhs != rhs;
}

inline bool is_after(const ContigRegion& lhs, const ContigRegion& rhs) noexcept
{
    return rhs.get_end() <= lhs.get_begin() && lhs != rhs;
}

inline bool are_adjacent(const ContigRegion& lhs, const ContigRegion& rhs) noexcept
{
    return lhs.get_begin() == rhs.get_end() || lhs.get_end() == rhs.get_begin();
}

inline ContigRegion::DifferenceType overlap_size(const ContigRegion& lhs, const ContigRegion& rhs) noexcept
{
    using DifferenceType = ContigRegion::DifferenceType;
    return static_cast<DifferenceType>(std::min(lhs.get_end(), rhs.get_end())) -
                    static_cast<DifferenceType>(std::max(lhs.get_begin(), rhs.get_begin()));
}

inline bool overlaps(const ContigRegion& lhs, const ContigRegion& rhs) noexcept
{
    const auto overlapped = overlap_size(lhs, rhs);
    return overlapped > 0 || (overlapped == 0 && (empty(lhs) || empty(rhs)));
}

inline bool contains(const ContigRegion& lhs, const ContigRegion& rhs) noexcept
{
    return lhs.get_begin() <= rhs.get_begin() && rhs.get_end() <= lhs.get_end();
}

inline ContigRegion::DifferenceType inner_distance(const ContigRegion& lhs,
                                                     const ContigRegion& rhs) noexcept
{
    if (overlaps(lhs, rhs)) return 0;
    
    return (is_before(lhs, rhs)) ? static_cast<ContigRegion::DifferenceType>(rhs.get_begin() - lhs.get_end())
                                : -inner_distance(rhs, lhs);
}

inline ContigRegion::DifferenceType outer_distance(const ContigRegion& lhs,
                                                     const ContigRegion& rhs) noexcept
{
    using DifferenceType = ContigRegion::DifferenceType;
    
    if (contains(lhs, rhs) || contains(rhs, lhs)) return 0;
    
    return static_cast<DifferenceType>(rhs.get_end()) - static_cast<DifferenceType>(lhs.get_begin());
}

inline ContigRegion shift(const ContigRegion& region, ContigRegion::DifferenceType n)
{
    using SizeType = ContigRegion::SizeType;
    
    if (n < 0 && region.get_begin() + n > region.get_begin()) {
        throw std::out_of_range {"shifted past contig start"};
    }
    
    return ContigRegion {
            static_cast<SizeType>(region.get_begin() + n), static_cast<SizeType>(region.get_end() + n)
    };
}

inline ContigRegion next_position(const ContigRegion& region)
{
    return ContigRegion {region.get_end(), region.get_end() + 1};
}

inline ContigRegion compress_lhs(const ContigRegion& region, ContigRegion::DifferenceType n)
{
    if (n < 0 && region.get_begin() + n > region.get_begin()) {
        throw std::out_of_range {"compressed past contig start"};
    }
    
    if (region.get_begin() + n > region.get_end()) {
        throw std::out_of_range {"compressed past region end"};
    }
    
    return ContigRegion {
        static_cast<ContigRegion::SizeType>(region.get_begin() + n),
        region.get_end()
    };
}

inline ContigRegion compress_rhs(const ContigRegion& region, ContigRegion::DifferenceType n)
{
    if (region.get_end() + n < region.get_begin()) {
        throw std::out_of_range {"compressed past region begin"};
    }
    
    return ContigRegion {
        region.get_begin(),
        static_cast<ContigRegion::SizeType>(region.get_end() + n)
    };
}

inline ContigRegion get_overlapped(const ContigRegion& lhs, const ContigRegion& rhs) noexcept
{
    if (!overlaps(lhs, rhs)) {
        throw std::runtime_error {"cannot get overlapped region between non overlapping regions"};
    }
    
    return ContigRegion {std::max(lhs.get_begin(), rhs.get_begin()), std::min(lhs.get_end(), rhs.get_end())};
}

inline ContigRegion get_encompassing(const ContigRegion& lhs, const ContigRegion& rhs) noexcept
{
    return ContigRegion {std::min(lhs.get_begin(), rhs.get_begin()), std::max(lhs.get_end(), rhs.get_end())};
}

inline ContigRegion get_intervening(const ContigRegion& lhs, const ContigRegion& rhs)
{
    if (begins_before(rhs, lhs) || overlaps(lhs, rhs)) {
        throw std::runtime_error {"cannot get intervening region between overlapping regions"};
    }
    
    return ContigRegion {lhs.get_end(), rhs.get_begin()};
}

inline ContigRegion::SizeType left_overhang_size(const ContigRegion& lhs, const ContigRegion& rhs) noexcept
{
    return (begins_before(lhs, rhs)) ? (rhs.get_begin() - lhs.get_begin()) : 0;
}

inline ContigRegion::SizeType right_overhang_size(const ContigRegion& lhs, const ContigRegion& rhs) noexcept
{
    return (ends_before(lhs, rhs)) ? 0 : (lhs.get_end() - rhs.get_end());
}

inline ContigRegion get_left_overhang(const ContigRegion& lhs, const ContigRegion& rhs) noexcept
{
    if (begins_before(rhs, lhs)) return ContigRegion {lhs.get_begin(), lhs.get_begin()};
    
    return ContigRegion {lhs.get_begin(), rhs.get_begin()};
}

inline ContigRegion get_right_overhang(const ContigRegion& lhs, const ContigRegion& rhs) noexcept
{
    if (ends_before(lhs, rhs)) return ContigRegion {lhs.get_end(), lhs.get_end()};
    
    return ContigRegion {rhs.get_end(), lhs.get_end()};
}

inline ContigRegion get_closed(const ContigRegion& lhs, const ContigRegion& rhs) noexcept
{
    return ContigRegion {lhs.get_begin(), rhs.get_end()};
}

inline ContigRegion get_head(const ContigRegion& region, ContigRegion::SizeType n = 0) noexcept
{
    const auto begin = region.get_begin();
    return ContigRegion {begin, std::min(begin + n, region.get_end())};
}

inline ContigRegion get_tail(const ContigRegion& region, ContigRegion::SizeType n = 0) noexcept
{
    const auto end = region.get_end();
    return ContigRegion {(end >= n) ? end - n : 0, end};
}

namespace std {
    template <> struct hash<ContigRegion>
    {
        size_t operator()(const ContigRegion& r) const
        {
            size_t seed {};
            boost::hash_combine(seed, r.get_begin());
            boost::hash_combine(seed, r.get_end());
            return seed;
        }
    };
}

inline std::string to_string(const ContigRegion& region)
{
    return std::to_string(region.get_begin()) + '-' + std::to_string(region.get_end());
}

inline std::ostream& operator<<(std::ostream& os, const ContigRegion& region)
{
    os << to_string(region);
    return os;
}

#endif /* defined(__Octopus__contig_region__) */
