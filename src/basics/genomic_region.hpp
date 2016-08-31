// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef genomic_region_hpp
#define genomic_region_hpp

#include <string>
#include <functional>
#include <stdexcept>
#include <ostream>
#include <cassert>

#include <boost/functional/hash.hpp>

#include <concepts/comparable.hpp>

#include "contig_region.hpp"

namespace octopus {

/**
    Represents a continuous region of a sequence in a genome. The contig
    name is the reference contig name (usually a chromosome), and the
    begin and end positions are zero-indexed half open - [begin,end) - indices.
 
    All comparison operations (<, ==, is_before, etc) throw exceptions if the arguements
    are not from the same contig.
*/
class GenomicRegion : public Comparable<GenomicRegion>
{
public:
    using ContigName = std::string;
    using Position   = ContigRegion::Position;
    using Size       = ContigRegion::Size;
    using Distance   = ContigRegion::Distance;
    
    GenomicRegion() = default;  // for use with containers
    
    template <typename T>
    explicit GenomicRegion(T&& contig_name, Position begin, Position end);
    
    template <typename T, typename R>
    explicit GenomicRegion(T&& contig_name, R&& contig_region);
    
    GenomicRegion(const GenomicRegion&)            = default;
    GenomicRegion& operator=(const GenomicRegion&) = default;
    GenomicRegion(GenomicRegion&&)                 = default;
    GenomicRegion& operator=(GenomicRegion&&)      = default;
    
    ~GenomicRegion() = default;
    
    const ContigName& contig_name() const noexcept;
    
    const ContigRegion& contig_region() const noexcept;
    
    Position begin() const noexcept;
    Position end() const noexcept;

private:
    ContigName contig_name_;
    ContigRegion contig_region_;
};

class RegionError : public std::logic_error
{
public:
    RegionError(GenomicRegion::ContigName first, GenomicRegion::ContigName second)
    :
    logic_error {"RegionError: comparings regions on different contigs"},
    first_ {first}, second_ {second} {}
    
    virtual ~RegionError() noexcept = default;
    
    virtual const char* what() const noexcept override
    {
        msg_ = std::string {logic_error::what()} + ": " + first_ + " & " + second_;
        return msg_.c_str();
    }
    
    auto first_contig() const { return first_; }
    auto second_contig() const { return second_; }
    
private:
    GenomicRegion::ContigName first_, second_;
    mutable std::string msg_;
};

// public member methods

template <typename T>
GenomicRegion::GenomicRegion(T&& contig_name, const Position begin, const Position end)
: contig_name_ {std::forward<T>(contig_name)}
, contig_region_ {begin, end}
{}

template <typename T, typename R>
GenomicRegion::GenomicRegion(T&& contig_name, R&& contig_region)
: contig_name_ {std::forward<T>(contig_name)}
, contig_region_ {std::forward<R>(contig_region)}
{}

inline const GenomicRegion::ContigName& GenomicRegion::contig_name() const noexcept
{
    return contig_name_;
}

inline const ContigRegion& GenomicRegion::contig_region() const noexcept
{
    return contig_region_;
}

inline GenomicRegion::Position GenomicRegion::begin() const noexcept
{
    return contig_region_.begin();
}

inline GenomicRegion::Position GenomicRegion::end() const noexcept
{
    return contig_region_.end();
}

// non-member methods

inline const GenomicRegion::ContigName& contig_name(const GenomicRegion& region) noexcept
{
    return region.contig_name();
}

inline GenomicRegion::Position mapped_begin(const GenomicRegion& region) noexcept
{
    return region.begin();
}

inline GenomicRegion::Position mapped_end(const GenomicRegion& region) noexcept
{
    return region.end();
}

inline std::string to_string(const GenomicRegion& region)
{
    return region.contig_name() + ':' + std::to_string(region.begin()) + '-'
                + std::to_string(region.end());
}

inline bool is_empty(const GenomicRegion& region) noexcept
{
    return is_empty(region.contig_region());
}

inline GenomicRegion::Size size(const GenomicRegion& region) noexcept
{
    return size(region.contig_region());
}

inline bool is_position(const GenomicRegion& region) noexcept
{
    return is_position(region.contig_region());
}

inline bool is_same_contig(const GenomicRegion& lhs, const GenomicRegion& rhs) noexcept
{
    return lhs.contig_name() == rhs.contig_name();
}

inline bool begins_equal(const GenomicRegion& lhs, const GenomicRegion& rhs)
{
    if (!is_same_contig(lhs, rhs)) throw RegionError {to_string(lhs), to_string(rhs)};
    return begins_equal(lhs.contig_region(), rhs.contig_region());
}

inline bool ends_equal(const GenomicRegion& lhs, const GenomicRegion& rhs)
{
    if (!is_same_contig(lhs, rhs)) throw RegionError {to_string(lhs), to_string(rhs)};
    return ends_equal(lhs.contig_region(), rhs.contig_region());
}

inline bool begins_before(const GenomicRegion& lhs, const GenomicRegion& rhs)
{
    if (!is_same_contig(lhs, rhs)) throw RegionError {to_string(lhs), to_string(rhs)};
    return begins_before(lhs.contig_region(), rhs.contig_region());
}

inline bool ends_before(const GenomicRegion& lhs, const GenomicRegion& rhs)
{
    if (!is_same_contig(lhs, rhs)) throw RegionError {to_string(lhs), to_string(rhs)};
    return ends_before(lhs.contig_region(), rhs.contig_region());
}

inline bool operator==(const GenomicRegion& lhs, const GenomicRegion& rhs) noexcept
{
    return is_same_contig(lhs, rhs) && lhs.contig_region() == rhs.contig_region();
}

inline bool operator<(const GenomicRegion& lhs, const GenomicRegion& rhs)
{
    if (!is_same_contig(lhs, rhs)) throw RegionError {to_string(lhs), to_string(rhs)};
    return lhs.contig_region() < rhs.contig_region();
}

inline bool is_before(const GenomicRegion& lhs, const GenomicRegion& rhs) noexcept
{
    if (!is_same_contig(lhs, rhs)) throw RegionError {to_string(lhs), to_string(rhs)};
    return is_before(lhs.contig_region(), rhs.contig_region());
}

inline bool is_after(const GenomicRegion& lhs, const GenomicRegion& rhs) noexcept
{
    if (!is_same_contig(lhs, rhs)) throw RegionError {to_string(lhs), to_string(rhs)};
    return is_after(lhs.contig_region(), rhs.contig_region());
}

inline bool are_adjacent(const GenomicRegion& lhs, const GenomicRegion& rhs) noexcept
{
    return is_same_contig(lhs, rhs) && are_adjacent(lhs.contig_region(), rhs.contig_region());
}

inline GenomicRegion::Distance overlap_size(const GenomicRegion& lhs, const GenomicRegion& rhs) noexcept
{
    return is_same_contig(lhs, rhs) ? overlap_size(lhs.contig_region(), rhs.contig_region()) : 0;
}

inline bool overlaps(const GenomicRegion& lhs, const GenomicRegion& rhs) noexcept
{
    return is_same_contig(lhs, rhs) && overlaps(lhs.contig_region(), rhs.contig_region());
}

inline bool contains(const GenomicRegion& lhs, const GenomicRegion& rhs) noexcept
{
    return is_same_contig(lhs, rhs) && contains(lhs.contig_region(), rhs.contig_region());
}

inline GenomicRegion::Distance inner_distance(const GenomicRegion& lhs, const GenomicRegion& rhs)
{
    if (!is_same_contig(lhs, rhs)) throw RegionError {to_string(lhs), to_string(rhs)};
    return inner_distance(lhs.contig_region(), rhs.contig_region());
}

inline GenomicRegion::Distance outer_distance(const GenomicRegion& lhs, const GenomicRegion& rhs)
{
    if (!is_same_contig(lhs, rhs)) throw RegionError {to_string(lhs), to_string(rhs)};
    return outer_distance(lhs.contig_region(), rhs.contig_region());
}

inline GenomicRegion shift(const GenomicRegion& region, GenomicRegion::Distance n)
{
    return GenomicRegion {region.contig_name(), shift(region.contig_region(), n)};
}

inline GenomicRegion next_position(const GenomicRegion& region)
{
    return GenomicRegion {region.contig_name(), next_position(region.contig_region())};
}

inline GenomicRegion expand_lhs(const GenomicRegion& region, const GenomicRegion::Distance n)
{
    return GenomicRegion {region.contig_name(), expand_lhs(region.contig_region(), n)};
}

inline GenomicRegion expand_rhs(const GenomicRegion& region, const GenomicRegion::Distance n)
{
    return GenomicRegion {region.contig_name(), expand_rhs(region.contig_region(), n)};
}

inline GenomicRegion expand(const GenomicRegion& region, const GenomicRegion::Distance n)
{
    return GenomicRegion {region.contig_name(), expand(region.contig_region(), n)};
}

inline GenomicRegion expand(const GenomicRegion& region, const GenomicRegion::Distance lhs,
                            const GenomicRegion::Distance rhs)
{
    return GenomicRegion {region.contig_name(), expand(region.contig_region(), lhs, rhs)};
}

inline GenomicRegion encompassing_region(const GenomicRegion& lhs, const GenomicRegion& rhs)
{
    if (!is_same_contig(lhs, rhs)) throw RegionError {to_string(lhs), to_string(rhs)};
    return GenomicRegion {lhs.contig_name(), encompassing_region(lhs.contig_region(), rhs.contig_region())};
}

inline GenomicRegion intervening_region(const GenomicRegion& lhs, const GenomicRegion& rhs)
{
    return GenomicRegion {lhs.contig_name(), intervening_region(lhs.contig_region(),  rhs.contig_region())};
}

inline GenomicRegion overlapped_region(const GenomicRegion& lhs, const GenomicRegion& rhs) noexcept
{
    return GenomicRegion {lhs.contig_name(), overlapped_region(lhs.contig_region(), rhs.contig_region())};
}

inline GenomicRegion::Size left_overhang_size(const GenomicRegion& lhs, const GenomicRegion& rhs) noexcept
{
    return (is_same_contig(lhs, rhs)) ? left_overhang_size(lhs.contig_region(), rhs.contig_region()) : 0;
}

inline GenomicRegion::Size right_overhang_size(const GenomicRegion& lhs, const GenomicRegion& rhs) noexcept
{
    return (is_same_contig(lhs, rhs)) ? right_overhang_size(lhs.contig_region(), rhs.contig_region()) : 0;
}

inline GenomicRegion left_overhang_region(const GenomicRegion& lhs, const GenomicRegion& rhs)
{
    if (!is_same_contig(lhs, rhs)) throw RegionError {to_string(lhs), to_string(rhs)};
    return GenomicRegion {lhs.contig_name(), left_overhang_region(lhs.contig_region(), rhs.contig_region())};
}

inline GenomicRegion right_overhang_region(const GenomicRegion& lhs, const GenomicRegion& rhs)
{
    if (!is_same_contig(lhs, rhs)) throw RegionError {to_string(lhs), to_string(rhs)};
    return GenomicRegion {lhs.contig_name(), right_overhang_region(lhs.contig_region(), rhs.contig_region())};
}

inline GenomicRegion closed_region(const GenomicRegion& lhs, const GenomicRegion& rhs)
{
    if (!is_same_contig(lhs, rhs)) throw RegionError {to_string(lhs), to_string(rhs)};
    return GenomicRegion {lhs.contig_name(), closed_region(lhs.contig_region(), rhs.contig_region())};
}

inline GenomicRegion head_region(const GenomicRegion& region, const GenomicRegion::Size n = 0)
{
    return GenomicRegion {region.contig_name(), head_region(region.contig_region(), n)};
}

inline GenomicRegion head_position(const GenomicRegion& region)
{
    return GenomicRegion {region.contig_name(), head_position(region.contig_region())};
}

inline GenomicRegion tail_region(const GenomicRegion& region, const GenomicRegion::Size n = 0)
{
    return GenomicRegion {region.contig_name(), tail_region(region.contig_region(), n)};
}

inline GenomicRegion tail_position(const GenomicRegion& region)
{
    return GenomicRegion {region.contig_name(), tail_position(region.contig_region())};
}

inline GenomicRegion::Distance begin_distance(const GenomicRegion& first, const GenomicRegion& second)
{
    if (!is_same_contig(first, second)) throw RegionError {to_string(first), to_string(second)};
    return begin_distance(first.contig_region(), second.contig_region());
}

inline GenomicRegion::Distance end_distance(const GenomicRegion& first, const GenomicRegion& second)
{
    if (!is_same_contig(first, second)) throw RegionError {to_string(first), to_string(second)};
    return end_distance(first.contig_region(), second.contig_region());
}

inline std::ostream& operator<<(std::ostream& os, const GenomicRegion& region)
{
    os << to_string(region);
    return os;
}
    
struct GenomicRegionHash
{
    std::size_t operator()(const GenomicRegion& region) const noexcept
    {
        using boost::hash_combine;
        std::size_t result {};
        hash_combine(result, std::hash<GenomicRegion::ContigName>()(region.contig_name()));
        hash_combine(result, std::hash<ContigRegion>()(region.contig_region()));
        return result;
    }
};
    
} // namespace octopus

namespace std {
    template <> struct hash<octopus::GenomicRegion>
    {
        size_t operator()(const octopus::GenomicRegion& region) const
        {
            return octopus::GenomicRegionHash()(region);
        }
    };
}

#endif
