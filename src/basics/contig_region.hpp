// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef contig_region_hpp
#define contig_region_hpp

#include <cstdint>
#include <algorithm>
#include <functional>
#include <cmath>
#include <stdexcept>
#include <ostream>

#include <boost/functional/hash.hpp>
#include <boost/optional.hpp>

#include "concepts/comparable.hpp"

namespace octopus {

/*
    Represents a region of continuous sequence.
    begin and end positions are zero-indexed half open indices [begin,end).
*/
class ContigRegion : public Comparable<ContigRegion>
{
public:
    using Position = std::uint_fast32_t;
    using Size     = Position;
    using Distance = std::int_fast64_t;
    
    class BadRegion;
    
    ContigRegion() = default;
    
    explicit ContigRegion(Position begin, Position end);
    
    ContigRegion(const ContigRegion&)            = default;
    ContigRegion& operator=(const ContigRegion&) = default;
    ContigRegion(ContigRegion&&)                 = default;
    ContigRegion& operator=(ContigRegion&&)      = default;
    
    ~ContigRegion() = default;
    
    Position begin() const noexcept { return begin_; }
    Position end() const noexcept { return end_; }
    
private:
    Position begin_, end_;
};

// BadRegion

class ContigRegion::BadRegion : public std::logic_error
{
public:
    using Position = ContigRegion::Position;
    
    BadRegion(Position begin, Position end) noexcept;
    
    virtual ~BadRegion() override = default;
    
    Position begin() const noexcept { return begin_; }
    Position end() const noexcept  { return end_; }
    
private:
    ContigRegion::Position begin_, end_;
};

// public member methods

inline ContigRegion::ContigRegion(const Position begin, const Position end)
: begin_ {begin}
, end_ {end}
{
    if (end < begin) throw BadRegion {begin, end};
}

inline ContigRegion::BadRegion::BadRegion(const Position begin, const Position end) noexcept
: std::logic_error {"BadRegion"}
, begin_ {begin}
, end_ {end}
{}

// non-member methods

inline ContigRegion::Position mapped_begin(const ContigRegion& region) noexcept
{
    return region.begin();
}

inline ContigRegion::Position mapped_end(const ContigRegion& region) noexcept
{
    return region.end();
}

inline bool is_empty(const ContigRegion& region) noexcept
{
    return region.begin() == region.end();
}

inline ContigRegion::Size size(const ContigRegion& region) noexcept
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

inline ContigRegion::Distance overlap_size(const ContigRegion& lhs, const ContigRegion& rhs) noexcept
{
    using D = ContigRegion::Distance;
    return static_cast<D>(std::min(lhs.end(), rhs.end()))
           - static_cast<D>(std::max(lhs.begin(), rhs.begin()));
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

inline ContigRegion::Distance inner_distance(const ContigRegion& lhs, const ContigRegion& rhs) noexcept
{
    if (overlaps(lhs, rhs)) return 0;
    return (is_before(lhs, rhs)) ? static_cast<ContigRegion::Distance>(rhs.begin() - lhs.end())
                                : -inner_distance(rhs, lhs);
}

inline ContigRegion::Distance outer_distance(const ContigRegion& lhs, const ContigRegion& rhs) noexcept
{
    using Distance = ContigRegion::Distance;
    if (contains(lhs, rhs) || contains(rhs, lhs)) return 0;
    return static_cast<Distance>(rhs.end()) - static_cast<Distance>(lhs.begin());
}

inline ContigRegion shift(const ContigRegion& region, ContigRegion::Distance n)
{
    using P = ContigRegion::Position;
    if (n < 0 && static_cast<P>(std::abs(n)) > region.begin()) {
        throw std::out_of_range {"ContigRegion: shifted past contig start"};
    }
    return ContigRegion {
        static_cast<P>(region.begin() + n),
        static_cast<P>(region.end() + n)
    };
}

inline ContigRegion next_position(const ContigRegion& region)
{
    return ContigRegion {region.end(), region.end() + 1};
}

inline ContigRegion expand_lhs(const ContigRegion& region, const ContigRegion::Distance n)
{
    using P = ContigRegion::Position;
    if (n > 0 && static_cast<P>(n) > region.begin()) {
        throw std::out_of_range {"ContigRegion: expanded past contig start"};
    }
    return ContigRegion {
        static_cast<P>(region.begin() - n),
        region.end()
    };
}

inline ContigRegion expand_rhs(const ContigRegion& region, const ContigRegion::Distance n)
{
    using P = ContigRegion::Position;
    if (n < 0 && static_cast<P>(std::abs(n)) > region.end()) {
        throw std::out_of_range {"ContigRegion: compressed past contig start"};
    }
    return ContigRegion {
        region.begin(),
        static_cast<P>(region.end() + n)
    };
}

inline ContigRegion expand(const ContigRegion& region, const ContigRegion::Distance n)
{
    using S = ContigRegion::Position;
    using D = ContigRegion::Distance;
    return ContigRegion {
        static_cast<S>(std::max(D {0}, static_cast<D>(region.begin()) - n)),
        static_cast<S>(region.end() + n)
    };
}

inline ContigRegion expand(const ContigRegion& region,
                           const ContigRegion::Distance lhs,
                           const ContigRegion::Distance rhs)
{
    return expand_lhs(expand_rhs(region, rhs), lhs);
}

inline boost::optional<ContigRegion> overlapped_region(const ContigRegion& lhs, const ContigRegion& rhs) noexcept
{
    if (!overlaps(lhs, rhs)) {
        return boost::none;
    }
    return ContigRegion {std::max(lhs.begin(), rhs.begin()), std::min(lhs.end(), rhs.end())};
}

inline ContigRegion encompassing_region(const ContigRegion& lhs, const ContigRegion& rhs) noexcept
{
    return ContigRegion {std::min(lhs.begin(), rhs.begin()), std::max(lhs.end(), rhs.end())};
}

inline boost::optional<ContigRegion> intervening_region(const ContigRegion& lhs, const ContigRegion& rhs) noexcept
{
    if (begins_before(rhs, lhs) || overlaps(lhs, rhs)) {
        return boost::none;
    }
    return ContigRegion {lhs.end(), rhs.begin()};
}

inline ContigRegion::Size left_overhang_size(const ContigRegion& lhs, const ContigRegion& rhs) noexcept
{
    return (begins_before(lhs, rhs)) ? (rhs.begin() - lhs.begin()) : 0;
}

inline ContigRegion::Size right_overhang_size(const ContigRegion& lhs, const ContigRegion& rhs) noexcept
{
    return (ends_before(lhs, rhs)) ? 0 : (lhs.end() - rhs.end());
}

inline ContigRegion left_overhang_region(const ContigRegion& lhs, const ContigRegion& rhs) noexcept
{
    if (begins_before(rhs, lhs)) {
        return ContigRegion {lhs.begin(), lhs.begin()};
    }
    return ContigRegion {lhs.begin(), rhs.begin()};
}

inline ContigRegion right_overhang_region(const ContigRegion& lhs, const ContigRegion& rhs) noexcept
{
    if (ends_before(lhs, rhs)) {
        return ContigRegion {lhs.end(), lhs.end()};
    }
    return ContigRegion {rhs.end(), lhs.end()};
}

inline ContigRegion closed_region(const ContigRegion& lhs, const ContigRegion& rhs) noexcept
{
    return ContigRegion {lhs.begin(), rhs.end()};
}

inline ContigRegion head_region(const ContigRegion& region, ContigRegion::Size n = 0) noexcept
{
    const auto begin = region.begin();
    return ContigRegion {begin, std::min(begin + n, region.end())};
}

inline ContigRegion head_position(const ContigRegion& region) noexcept
{
    return head_region(region, 1);
}

inline ContigRegion tail_region(const ContigRegion& region, const ContigRegion::Size n = 0) noexcept
{
    const auto end = region.end();
    return ContigRegion {(end >= n) ? end - n : 0, end};
}

inline ContigRegion tail_position(const ContigRegion& region) noexcept
{
    return tail_region(region, 1);
}

inline ContigRegion::Distance begin_distance(const ContigRegion& first, const ContigRegion& second) noexcept
{
    return static_cast<ContigRegion::Distance>(second.begin()) - first.begin();
}

inline ContigRegion::Distance end_distance(const ContigRegion& first, const ContigRegion& second) noexcept
{
    return static_cast<ContigRegion::Distance>(second.end()) - first.end();
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

struct ContigRegionHash
{
    std::size_t operator()(const ContigRegion& region) const noexcept
    {
        using boost::hash_combine;
        std::size_t seed {};
        hash_combine(seed, region.begin());
        hash_combine(seed, region.end());
        return seed;
    }
};

} // namespace octopus

namespace std {
    template <> struct hash<octopus::ContigRegion>
    {
        size_t operator()(const octopus::ContigRegion& region) const noexcept
        {
            return octopus::ContigRegionHash()(region);
        }
    };
}

#endif
