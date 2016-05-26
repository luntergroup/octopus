//
//  genomic_region.hpp
//  Octopus
//
//  Created by Daniel Cooke on 09/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__genomic_region__
#define __Octopus__genomic_region__

#include <string>
#include <ostream>
#include <stdexcept>
#include <vector>
#include <cassert>

#include <boost/functional/hash.hpp>

#include "contig_region.hpp"
#include "comparable.hpp"

/**
    Represents a continuous region of a sequence in a genome. The sequence
    name is the reference sequence name (usually a chromosome), and the
    begin and end positions are zero-indexed half open - [begin,end) - indices.
 
    All comparison operations (<, ==, is_before, etc) throw exceptions if the arguements
    are not from the same contig.
*/
class GenomicRegion : public Comparable<GenomicRegion>
{
public:
    using ContigNameType = std::string;
    using SizeType       = ContigRegion::SizeType;
    using DifferenceType = ContigRegion::DifferenceType;
    
    GenomicRegion() = default;  // for use with containers
    
    template <typename T> explicit GenomicRegion(T&& contig_name, SizeType begin, SizeType end);
    
    template <typename T, typename R> explicit GenomicRegion(T&& contig_name, R&& contig_region);
    
    ~GenomicRegion() = default;
    
    GenomicRegion(const GenomicRegion&)            = default;
    GenomicRegion& operator=(const GenomicRegion&) = default;
    GenomicRegion(GenomicRegion&&)                 = default;
    GenomicRegion& operator=(GenomicRegion&&)      = default;
    
    const ContigNameType& contig_name() const noexcept;
    const ContigRegion& contig_region() const noexcept;
    SizeType begin() const noexcept;
    SizeType end() const noexcept;

private:
    ContigNameType contig_name_;
    ContigRegion contig_region_;
};

class RegionError : public std::runtime_error
{
public:
    RegionError(GenomicRegion::ContigNameType first, GenomicRegion::ContigNameType second)
    :
    runtime_error {"cannot compare regions on different contigs"},
    first_ {first}, second_ {second} {}
    
    const char* what() const noexcept
    {
        return (std::string{runtime_error::what()} + ": " + first_ + " & " + second_).c_str();
    }
    
private:
    GenomicRegion::ContigNameType first_, second_;
};

// public member methods

template <typename T>
GenomicRegion::GenomicRegion(T&& contig_name, const SizeType begin, const SizeType end)
:
contig_name_ {std::forward<T>(contig_name)},
contig_region_ {begin, end}
{}

template <typename T, typename R>
GenomicRegion::GenomicRegion(T&& contig_name, R&& contig_region)
:
contig_name_ {std::forward<T>(contig_name)},
contig_region_ {std::forward<R>(contig_region)}
{}

inline const GenomicRegion::ContigNameType& GenomicRegion::contig_name() const noexcept
{
    return contig_name_;
}

inline const ContigRegion& GenomicRegion::contig_region() const noexcept
{
    return contig_region_;
}

inline GenomicRegion::SizeType GenomicRegion::begin() const noexcept
{
    return contig_region_.begin();
}

inline GenomicRegion::SizeType GenomicRegion::end() const noexcept
{
    return contig_region_.end();
}

// non-member methods

inline const GenomicRegion::ContigNameType& contig_name(const GenomicRegion& region) noexcept
{
    return region.contig_name();
}

inline GenomicRegion::SizeType mapped_begin(const GenomicRegion& region) noexcept
{
    return region.begin();
}

inline GenomicRegion::SizeType mapped_end(const GenomicRegion& region) noexcept
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

inline GenomicRegion::SizeType size(const GenomicRegion& region) noexcept
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

inline GenomicRegion::DifferenceType overlap_size(const GenomicRegion& lhs, const GenomicRegion& rhs) noexcept
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

inline GenomicRegion::DifferenceType inner_distance(const GenomicRegion& lhs, const GenomicRegion& rhs)
{
    if (!is_same_contig(lhs, rhs)) throw RegionError {to_string(lhs), to_string(rhs)};
    return inner_distance(lhs.contig_region(), rhs.contig_region());
}

inline GenomicRegion::DifferenceType outer_distance(const GenomicRegion& lhs, const GenomicRegion& rhs)
{
    if (!is_same_contig(lhs, rhs)) throw RegionError {to_string(lhs), to_string(rhs)};
    return outer_distance(lhs.contig_region(), rhs.contig_region());
}

inline GenomicRegion shift(const GenomicRegion& region, GenomicRegion::DifferenceType n)
{
    return GenomicRegion {region.contig_name(), shift(region.contig_region(), n)};
}

inline GenomicRegion next_position(const GenomicRegion& region)
{
    return GenomicRegion {region.contig_name(), next_position(region.contig_region())};
}

inline GenomicRegion expand_lhs(const GenomicRegion& region, const GenomicRegion::DifferenceType n)
{
    return GenomicRegion {region.contig_name(), expand_lhs(region.contig_region(), n)};
}

inline GenomicRegion expand_rhs(const GenomicRegion& region, const GenomicRegion::DifferenceType n)
{
    return GenomicRegion {region.contig_name(), expand_rhs(region.contig_region(), n)};
}

inline GenomicRegion expand(const GenomicRegion& region, const GenomicRegion::DifferenceType n)
{
    return GenomicRegion {region.contig_name(), expand(region.contig_region(), n)};
}

inline GenomicRegion encompassing_region(const GenomicRegion& lhs, const GenomicRegion& rhs)
{
    if (!is_same_contig(lhs, rhs)) throw RegionError {to_string(lhs), to_string(rhs)};
    return GenomicRegion {lhs.contig_name(), encompassing_region(lhs.contig_region(),
                                                                     rhs.contig_region())};
}

inline GenomicRegion intervening_region(const GenomicRegion& lhs, const GenomicRegion& rhs)
{
    return GenomicRegion {lhs.contig_name(), intervening_region(lhs.contig_region(),
                                                                    rhs.contig_region())};
}

inline GenomicRegion overlapped_region(const GenomicRegion& lhs, const GenomicRegion& rhs) noexcept
{
    return GenomicRegion {lhs.contig_name(), overlapped_region(lhs.contig_region(),
                                                                   rhs.contig_region())};
}

inline GenomicRegion::SizeType left_overhang_size(const GenomicRegion& lhs, const GenomicRegion& rhs) noexcept
{
    return (is_same_contig(lhs, rhs)) ? left_overhang_size(lhs.contig_region(), rhs.contig_region()) : 0;
}

inline GenomicRegion::SizeType right_overhang_size(const GenomicRegion& lhs, const GenomicRegion& rhs) noexcept
{
    return (is_same_contig(lhs, rhs)) ? right_overhang_size(lhs.contig_region(), rhs.contig_region()) : 0;
}

inline GenomicRegion left_overhang_region(const GenomicRegion& lhs, const GenomicRegion& rhs)
{
    if (!is_same_contig(lhs, rhs)) throw RegionError {to_string(lhs), to_string(rhs)};
    return GenomicRegion {lhs.contig_name(), left_overhang_region(lhs.contig_region(),
                                                                      rhs.contig_region())};
}

inline GenomicRegion right_overhang_region(const GenomicRegion& lhs, const GenomicRegion& rhs)
{
    if (!is_same_contig(lhs, rhs)) throw RegionError {to_string(lhs), to_string(rhs)};
    return GenomicRegion {lhs.contig_name(), right_overhang_region(lhs.contig_region(),
                                                                       rhs.contig_region())};
}

inline GenomicRegion closed_region(const GenomicRegion& lhs, const GenomicRegion& rhs)
{
    if (!is_same_contig(lhs, rhs)) throw RegionError {to_string(lhs), to_string(rhs)};
    return GenomicRegion {lhs.contig_name(), closed_region(lhs.contig_region(),
                                                               rhs.contig_region())};
}

inline GenomicRegion head_region(const GenomicRegion& region, const GenomicRegion::SizeType n = 0)
{
    return GenomicRegion {region.contig_name(), head_region(region.contig_region(), n)};
}

inline GenomicRegion head_position(const GenomicRegion& region)
{
    return GenomicRegion {region.contig_name(), head_position(region.contig_region())};
}

inline GenomicRegion tail_region(const GenomicRegion& region, const GenomicRegion::SizeType n = 0)
{
    return GenomicRegion {region.contig_name(), tail_region(region.contig_region(), n)};
}

inline GenomicRegion tail_position(const GenomicRegion& region)
{
    return GenomicRegion {region.contig_name(), tail_position(region.contig_region())};
}

inline GenomicRegion::DifferenceType begin_distance(const GenomicRegion& lhs, const GenomicRegion& rhs)
{
    if (!is_same_contig(lhs, rhs)) throw RegionError {to_string(lhs), to_string(rhs)};
    return begin_distance(lhs.contig_region(), rhs.contig_region());
}

inline GenomicRegion::DifferenceType end_distance(const GenomicRegion& lhs, const GenomicRegion& rhs)
{
    if (!is_same_contig(lhs, rhs)) throw RegionError {to_string(lhs), to_string(rhs)};
    return end_distance(lhs.contig_region(), rhs.contig_region());
}

namespace std {
    template <> struct hash<GenomicRegion>
    {
        size_t operator()(const GenomicRegion& r) const
        {
            size_t result {};
            boost::hash_combine(result, hash<GenomicRegion::ContigNameType>()(r.contig_name()));
            boost::hash_combine(result, hash<ContigRegion>()(r.contig_region()));
            return result;
        }
    };
}

inline std::ostream& operator<<(std::ostream& os, const GenomicRegion& region)
{
    os << to_string(region);
    return os;
}

#endif /* defined(__Octopus__genomic_region__) */
