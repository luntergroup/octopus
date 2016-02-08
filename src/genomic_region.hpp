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
    
    GenomicRegion() = default;
    template <typename T> explicit GenomicRegion(T&& contig_name, SizeType begin, SizeType end);
    template <typename T, typename R> explicit GenomicRegion(T&& contig_name, R&& contig_region);
    ~GenomicRegion() = default;
    
    GenomicRegion(const GenomicRegion&)            = default;
    GenomicRegion& operator=(const GenomicRegion&) = default;
    GenomicRegion(GenomicRegion&&)                 = default;
    GenomicRegion& operator=(GenomicRegion&&)      = default;
    
    const ContigNameType& get_contig_name() const noexcept;
    const ContigRegion& get_contig_region() const noexcept;
    SizeType get_begin() const noexcept;
    SizeType get_end() const noexcept;

private:
    ContigNameType contig_name_;
    ContigRegion contig_region_;
};

class RegionError : std::runtime_error
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
GenomicRegion::GenomicRegion(T&& contig_name, SizeType begin, SizeType end)
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

inline const GenomicRegion::ContigNameType& GenomicRegion::get_contig_name() const noexcept
{
    return contig_name_;
}

inline const ContigRegion& GenomicRegion::get_contig_region() const noexcept
{
    return contig_region_;
}

inline GenomicRegion::SizeType GenomicRegion::get_begin() const noexcept
{
    return contig_region_.get_begin();
}

inline GenomicRegion::SizeType GenomicRegion::get_end() const noexcept
{
    return contig_region_.get_end();
}

// non-member methods

inline const GenomicRegion::ContigNameType& get_contig_name(const GenomicRegion& region) noexcept
{
    return region.get_contig_name();
}

inline GenomicRegion::SizeType get_begin(const GenomicRegion& region) noexcept
{
    return region.get_begin();
}

inline GenomicRegion::SizeType get_end(const GenomicRegion& region) noexcept
{
    return region.get_end();
}

inline std::string to_string(const GenomicRegion& region)
{
    return region.get_contig_name() + ':' + std::to_string(region.get_begin()) + '-'
    + std::to_string(region.get_end());
}

inline bool is_empty(const GenomicRegion& region) noexcept
{
    return is_empty(region.get_contig_region());
}

inline GenomicRegion::SizeType region_size(const GenomicRegion& region) noexcept
{
    return region_size(region.get_contig_region());
}

inline bool is_same_contig(const GenomicRegion& lhs, const GenomicRegion& rhs) noexcept
{
    return lhs.get_contig_name() == rhs.get_contig_name();
}

inline bool begins_equal(const GenomicRegion& lhs, const GenomicRegion& rhs)
{
    if (!is_same_contig(lhs, rhs)) throw RegionError {to_string(lhs), to_string(rhs)};
    return begins_equal(lhs.get_contig_region(), rhs.get_contig_region());
}

inline bool ends_equal(const GenomicRegion& lhs, const GenomicRegion& rhs)
{
    if (!is_same_contig(lhs, rhs)) throw RegionError {to_string(lhs), to_string(rhs)};
    return ends_equal(lhs.get_contig_region(), rhs.get_contig_region());
}

inline bool begins_before(const GenomicRegion& lhs, const GenomicRegion& rhs)
{
    if (!is_same_contig(lhs, rhs)) throw RegionError {to_string(lhs), to_string(rhs)};
    return begins_before(lhs.get_contig_region(), rhs.get_contig_region());
}

inline bool ends_before(const GenomicRegion& lhs, const GenomicRegion& rhs)
{
    if (!is_same_contig(lhs, rhs)) throw RegionError {to_string(lhs), to_string(rhs)};
    return ends_before(lhs.get_contig_region(), rhs.get_contig_region());
}

inline bool operator==(const GenomicRegion& lhs, const GenomicRegion& rhs) noexcept
{
    return is_same_contig(lhs, rhs) && lhs.get_contig_region() == rhs.get_contig_region();
}

inline bool operator<(const GenomicRegion& lhs, const GenomicRegion& rhs)
{
    if (!is_same_contig(lhs, rhs)) throw RegionError {to_string(lhs), to_string(rhs)};
    return lhs.get_contig_region() < rhs.get_contig_region();
}

inline bool is_before(const GenomicRegion& lhs, const GenomicRegion& rhs) noexcept
{
    if (!is_same_contig(lhs, rhs)) throw RegionError {to_string(lhs), to_string(rhs)};
    return is_before(lhs.get_contig_region(), rhs.get_contig_region());
}

inline bool is_after(const GenomicRegion& lhs, const GenomicRegion& rhs) noexcept
{
    if (!is_same_contig(lhs, rhs)) throw RegionError {to_string(lhs), to_string(rhs)};
    return is_after(lhs.get_contig_region(), rhs.get_contig_region());
}

inline bool are_adjacent(const GenomicRegion& lhs, const GenomicRegion& rhs) noexcept
{
    return is_same_contig(lhs, rhs) && are_adjacent(lhs.get_contig_region(), rhs.get_contig_region());
}

inline GenomicRegion::DifferenceType overlap_size(const GenomicRegion& lhs, const GenomicRegion& rhs) noexcept
{
    return is_same_contig(lhs, rhs) ? overlap_size(lhs.get_contig_region(), rhs.get_contig_region()) : 0;
}

inline bool overlaps(const GenomicRegion& lhs, const GenomicRegion& rhs) noexcept
{
    return is_same_contig(lhs, rhs) && overlaps(lhs.get_contig_region(), rhs.get_contig_region());
}

inline bool contains(const GenomicRegion& lhs, const GenomicRegion& rhs) noexcept
{
    return is_same_contig(lhs, rhs) && contains(lhs.get_contig_region(), rhs.get_contig_region());
}

inline GenomicRegion::DifferenceType inner_distance(const GenomicRegion& lhs, const GenomicRegion& rhs)
{
    if (!is_same_contig(lhs, rhs)) throw RegionError {to_string(lhs), to_string(rhs)};
    return inner_distance(lhs.get_contig_region(), rhs.get_contig_region());
}

inline GenomicRegion::DifferenceType outer_distance(const GenomicRegion& lhs, const GenomicRegion& rhs)
{
    if (!is_same_contig(lhs, rhs)) throw RegionError {to_string(lhs), to_string(rhs)};
    return outer_distance(lhs.get_contig_region(), rhs.get_contig_region());
}

inline GenomicRegion shift(const GenomicRegion& region, GenomicRegion::DifferenceType n)
{
    return GenomicRegion {region.get_contig_name(), shift(region.get_contig_region(), n)};
}

inline GenomicRegion next_position(const GenomicRegion& region)
{
    return GenomicRegion {region.get_contig_name(), next_position(region.get_contig_region())};
}

inline GenomicRegion compress_lhs(const GenomicRegion& region, GenomicRegion::DifferenceType n)
{
    return GenomicRegion {region.get_contig_name(), compress_lhs(region.get_contig_region(), n)};
}

inline GenomicRegion compress_rhs(const GenomicRegion& region, GenomicRegion::DifferenceType n)
{
    return GenomicRegion {region.get_contig_name(), compress_rhs(region.get_contig_region(), n)};
}

inline GenomicRegion encompassing_region(const GenomicRegion& lhs, const GenomicRegion& rhs)
{
    if (!is_same_contig(lhs, rhs)) throw RegionError {to_string(lhs), to_string(rhs)};
    return GenomicRegion {lhs.get_contig_name(), encompassing_region(lhs.get_contig_region(),
                                                                     rhs.get_contig_region())};
}

inline GenomicRegion intervening_region(const GenomicRegion& lhs, const GenomicRegion& rhs)
{
    return GenomicRegion {lhs.get_contig_name(), intervening_region(lhs.get_contig_region(),
                                                                    rhs.get_contig_region())};
}

inline GenomicRegion overlapped_region(const GenomicRegion& lhs, const GenomicRegion& rhs) noexcept
{
    return GenomicRegion {lhs.get_contig_name(), overlapped_region(lhs.get_contig_region(),
                                                                   rhs.get_contig_region())};
}

inline GenomicRegion::SizeType left_overhang_size(const GenomicRegion& lhs, const GenomicRegion& rhs) noexcept
{
    return (is_same_contig(lhs, rhs)) ? left_overhang_size(lhs.get_contig_region(), rhs.get_contig_region()) : 0;
}

inline GenomicRegion::SizeType right_overhang_size(const GenomicRegion& lhs, const GenomicRegion& rhs) noexcept
{
    return (is_same_contig(lhs, rhs)) ? right_overhang_size(lhs.get_contig_region(), rhs.get_contig_region()) : 0;
}

inline GenomicRegion left_overhang_region(const GenomicRegion& lhs, const GenomicRegion& rhs)
{
    if (!is_same_contig(lhs, rhs)) throw RegionError {to_string(lhs), to_string(rhs)};
    return GenomicRegion {lhs.get_contig_name(), left_overhang_region(lhs.get_contig_region(),
                                                                      rhs.get_contig_region())};
}

inline GenomicRegion right_overhang_region(const GenomicRegion& lhs, const GenomicRegion& rhs)
{
    if (!is_same_contig(lhs, rhs)) throw RegionError {to_string(lhs), to_string(rhs)};
    return GenomicRegion {lhs.get_contig_name(), right_overhang_region(lhs.get_contig_region(),
                                                                       rhs.get_contig_region())};
}

inline GenomicRegion closed_region(const GenomicRegion& lhs, const GenomicRegion& rhs)
{
    if (!is_same_contig(lhs, rhs)) throw RegionError {to_string(lhs), to_string(rhs)};
    return GenomicRegion {lhs.get_contig_name(), closed_region(lhs.get_contig_region(),
                                                               rhs.get_contig_region())};
}

inline GenomicRegion head_region(const GenomicRegion& region, GenomicRegion::SizeType n = 0)
{
    return GenomicRegion {region.get_contig_name(), head_region(region.get_contig_region(), n)};
}

inline GenomicRegion head_position(const GenomicRegion& region)
{
    return GenomicRegion {region.get_contig_name(), head_position(region.get_contig_region())};
}

inline GenomicRegion tail_region(const GenomicRegion& region, GenomicRegion::SizeType n = 0)
{
    return GenomicRegion {region.get_contig_name(), tail_region(region.get_contig_region(), n)};
}

inline GenomicRegion tail_position(const GenomicRegion& region)
{
    return GenomicRegion {region.get_contig_name(), tail_position(region.get_contig_region())};
}

inline GenomicRegion::DifferenceType begin_distance(const GenomicRegion& lhs, const GenomicRegion& rhs)
{
    if (!is_same_contig(lhs, rhs)) throw RegionError {to_string(lhs), to_string(rhs)};
    return begin_distance(lhs.get_contig_region(), rhs.get_contig_region());
}

inline GenomicRegion::DifferenceType end_distance(const GenomicRegion& lhs, const GenomicRegion& rhs)
{
    if (!is_same_contig(lhs, rhs)) throw RegionError {to_string(lhs), to_string(rhs)};
    return end_distance(lhs.get_contig_region(), rhs.get_contig_region());
}

namespace std {
    template <> struct hash<GenomicRegion>
    {
        size_t operator()(const GenomicRegion& r) const
        {
            size_t result {};
            boost::hash_combine(result, hash<GenomicRegion::ContigNameType>()(r.get_contig_name()));
            boost::hash_combine(result, hash<ContigRegion>()(r.get_contig_region()));
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
