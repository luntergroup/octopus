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
    
    SizeType begin() const noexcept;
    SizeType end() const noexcept;
    
private:
    SizeType begin_, end_;
};

// public member methods

inline ContigRegion::ContigRegion(const SizeType begin, const SizeType end)
:
begin_ {begin},
end_ {end}
{
    if (end < begin) throw std::runtime_error {"ContigRegion constructed with end < begin"};
}

inline ContigRegion::SizeType ContigRegion::begin() const noexcept
{
    return begin_;
}

inline ContigRegion::SizeType ContigRegion::end() const noexcept
{
    return end_;
}

// non-member methods

inline ContigRegion::SizeType mapped_begin(const ContigRegion& region) noexcept
{
    return region.begin();
}

inline ContigRegion::SizeType mapped_end(const ContigRegion& region) noexcept
{
    return region.end();
}

inline bool is_empty(const ContigRegion& region) noexcept
{
    return region.begin() == region.end();
}

inline ContigRegion::SizeType size(const ContigRegion& region) noexcept
{
    return region.end() - region.begin();
}

inline bool is_position(const ContigRegion& region) noexcept
{
    return size(region) == 1;
}

inline bool begins_equal(const ContigRegion& lhs, const ContigRegion& rhs) noexcept
{
    return lhs.begin() == rhs.begin();
}

inline bool ends_equal(const ContigRegion& lhs, const ContigRegion& rhs) noexcept
{
    return lhs.end() == rhs.end();
}

inline bool begins_before(const ContigRegion& lhs, const ContigRegion& rhs) noexcept
{
    return lhs.begin() < rhs.begin();
}

inline bool ends_before(const ContigRegion& lhs, const ContigRegion& rhs) noexcept
{
    return lhs.end() < rhs.end();
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
    return lhs.end() <= rhs.begin() && lhs != rhs;
}

inline bool is_after(const ContigRegion& lhs, const ContigRegion& rhs) noexcept
{
    return rhs.end() <= lhs.begin() && lhs != rhs;
}

inline bool are_adjacent(const ContigRegion& lhs, const ContigRegion& rhs) noexcept
{
    return lhs.begin() == rhs.end() || lhs.end() == rhs.begin();
}

inline ContigRegion::DifferenceType overlap_size(const ContigRegion& lhs, const ContigRegion& rhs) noexcept
{
    using DifferenceType = ContigRegion::DifferenceType;
    return static_cast<DifferenceType>(std::min(lhs.end(), rhs.end())) -
                    static_cast<DifferenceType>(std::max(lhs.begin(), rhs.begin()));
}

inline bool overlaps(const ContigRegion& lhs, const ContigRegion& rhs) noexcept
{
    const auto overlapped = overlap_size(lhs, rhs);
    return overlapped > 0 || (overlapped == 0 && (is_empty(lhs) || is_empty(rhs)));
}

inline bool contains(const ContigRegion& lhs, const ContigRegion& rhs) noexcept
{
    return lhs.begin() <= rhs.begin() && rhs.end() <= lhs.end();
}

inline ContigRegion::DifferenceType inner_distance(const ContigRegion& lhs,
                                                     const ContigRegion& rhs) noexcept
{
    if (overlaps(lhs, rhs)) return 0;
    
    return (is_before(lhs, rhs)) ? static_cast<ContigRegion::DifferenceType>(rhs.begin() - lhs.end())
                                : -inner_distance(rhs, lhs);
}

inline ContigRegion::DifferenceType outer_distance(const ContigRegion& lhs,
                                                     const ContigRegion& rhs) noexcept
{
    using DifferenceType = ContigRegion::DifferenceType;
    
    if (contains(lhs, rhs) || contains(rhs, lhs)) return 0;
    
    return static_cast<DifferenceType>(rhs.end()) - static_cast<DifferenceType>(lhs.begin());
}

inline ContigRegion shift(const ContigRegion& region, ContigRegion::DifferenceType n)
{
    using SizeType = ContigRegion::SizeType;
    
    if (n < 0 && region.begin() + n > region.begin()) {
        throw std::out_of_range {"shifted past contig start"};
    }
    
    return ContigRegion {
            static_cast<SizeType>(region.begin() + n), static_cast<SizeType>(region.end() + n)
    };
}

inline ContigRegion next_position(const ContigRegion& region)
{
    return ContigRegion {region.end(), region.end() + 1};
}

inline ContigRegion expand_lhs(const ContigRegion& region, const ContigRegion::DifferenceType n)
{
    if (n < 0 && region.begin() + n > region.begin()) {
        throw std::out_of_range {"compressed past contig start"};
    }
    
    if (region.begin() - n > region.end()) {
        throw std::out_of_range {"compressed past region end"};
    }
    
    return ContigRegion {
        static_cast<ContigRegion::SizeType>(region.begin() - n),
        region.end()
    };
}

inline ContigRegion expand_rhs(const ContigRegion& region, const ContigRegion::DifferenceType n)
{
    if (region.end() + n < region.begin()) {
        throw std::out_of_range {"compressed past region begin"};
    }
    
    return ContigRegion {
        region.begin(),
        static_cast<ContigRegion::SizeType>(region.end() + n)
    };
}

inline ContigRegion expand(const ContigRegion& region, const ContigRegion::DifferenceType n)
{
    using S = ContigRegion::SizeType;
    using D = ContigRegion::DifferenceType;
    return ContigRegion {
        static_cast<S>(std::max(D {0}, static_cast<D>(region.begin()) - n)),
        static_cast<S>(region.end() + n)
    };
}

inline ContigRegion overlapped_region(const ContigRegion& lhs, const ContigRegion& rhs) noexcept
{
    if (!overlaps(lhs, rhs)) {
        throw std::runtime_error {"cannot calculate overlapped region between non overlapping regions"};
    }
    
    return ContigRegion {std::max(lhs.begin(), rhs.begin()), std::min(lhs.end(), rhs.end())};
}

inline ContigRegion encompassing_region(const ContigRegion& lhs, const ContigRegion& rhs) noexcept
{
    return ContigRegion {std::min(lhs.begin(), rhs.begin()), std::max(lhs.end(), rhs.end())};
}

inline ContigRegion intervening_region(const ContigRegion& lhs, const ContigRegion& rhs)
{
    if (begins_before(rhs, lhs) || overlaps(lhs, rhs)) {
        throw std::runtime_error {"cannot calculate intervening region between overlapping regions"};
    }
    
    return ContigRegion {lhs.end(), rhs.begin()};
}

inline ContigRegion::SizeType left_overhang_size(const ContigRegion& lhs, const ContigRegion& rhs) noexcept
{
    return (begins_before(lhs, rhs)) ? (rhs.begin() - lhs.begin()) : 0;
}

inline ContigRegion::SizeType right_overhang_size(const ContigRegion& lhs, const ContigRegion& rhs) noexcept
{
    return (ends_before(lhs, rhs)) ? 0 : (lhs.end() - rhs.end());
}

inline ContigRegion left_overhang_region(const ContigRegion& lhs, const ContigRegion& rhs) noexcept
{
    if (begins_before(rhs, lhs)) return ContigRegion {lhs.begin(), lhs.begin()};
    
    return ContigRegion {lhs.begin(), rhs.begin()};
}

inline ContigRegion right_overhang_region(const ContigRegion& lhs, const ContigRegion& rhs) noexcept
{
    if (ends_before(lhs, rhs)) return ContigRegion {lhs.end(), lhs.end()};
    
    return ContigRegion {rhs.end(), lhs.end()};
}

inline ContigRegion closed_region(const ContigRegion& lhs, const ContigRegion& rhs) noexcept
{
    return ContigRegion {lhs.begin(), rhs.end()};
}

inline ContigRegion head_region(const ContigRegion& region, ContigRegion::SizeType n = 0) noexcept
{
    const auto begin = region.begin();
    return ContigRegion {begin, std::min(begin + n, region.end())};
}

inline ContigRegion head_position(const ContigRegion& region) noexcept
{
    return head_region(region, 1);
}

inline ContigRegion tail_region(const ContigRegion& region, const ContigRegion::SizeType n = 0) noexcept
{
    const auto end = region.end();
    return ContigRegion {(end >= n) ? end - n : 0, end};
}

inline ContigRegion tail_position(const ContigRegion& region) noexcept
{
    return tail_region(region, 1);
}

inline ContigRegion::DifferenceType begin_distance(const ContigRegion& lhs, const ContigRegion& rhs) noexcept
{
    return static_cast<ContigRegion::DifferenceType>(lhs.begin()) - rhs.begin();
}

inline ContigRegion::DifferenceType end_distance(const ContigRegion& lhs, const ContigRegion& rhs) noexcept
{
    return static_cast<ContigRegion::DifferenceType>(lhs.end()) - rhs.end();
}

namespace std {
    template <> struct hash<ContigRegion>
    {
        size_t operator()(const ContigRegion& r) const
        {
            size_t seed {};
            boost::hash_combine(seed, r.begin());
            boost::hash_combine(seed, r.end());
            return seed;
        }
    };
}

inline std::string to_string(const ContigRegion& region)
{
    return std::to_string(region.begin()) + '-' + std::to_string(region.end());
}

inline std::ostream& operator<<(std::ostream& os, const ContigRegion& region)
{
    os << to_string(region);
    return os;
}

#endif /* defined(__Octopus__contig_region__) */
