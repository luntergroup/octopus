//
//  genomic_region.h
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
#include <boost/functional/hash.hpp> // boost::hash_combine

#include "contig_region.h"
#include "comparable.h"

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
    using StringType     = std::string;
    using SizeType       = ContigRegion::SizeType;
    using DifferenceType = ContigRegion::DifferenceType;
    
    GenomicRegion() = default;
    template <typename T>
    explicit GenomicRegion(T&& contig_name, SizeType begin, SizeType end);
    template <typename T, typename R>
    explicit GenomicRegion(T&& contig_name, R&& contig_region);
    
    GenomicRegion(const GenomicRegion&)            = default;
    GenomicRegion& operator=(const GenomicRegion&) = default;
    GenomicRegion(GenomicRegion&&)                 = default;
    GenomicRegion& operator=(GenomicRegion&&)      = default;
    
    const StringType& get_contig_name() const noexcept;
    const ContigRegion& get_contig_region() const noexcept;
    SizeType get_begin() const noexcept;
    SizeType get_end() const noexcept;

private:
    StringType contig_name_;
    ContigRegion contig_region_;
};

class RegionError : std::runtime_error {
public:
    RegionError(GenomicRegion::StringType first, GenomicRegion::StringType second)
    :
    runtime_error {"Cannot compare regions on different contigs"},
    first_ {first}, second_ {second} {}
    
    const char* what() const noexcept
    {
        return (std::string{runtime_error::what()} + ": " + first_ + " & " + second_).c_str();
    }
    
private:
    GenomicRegion::StringType first_, second_;
};

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

inline const GenomicRegion::StringType& GenomicRegion::get_contig_name() const noexcept
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

inline const GenomicRegion::StringType& get_contig_name(const GenomicRegion& a_region) noexcept
{
    return a_region.get_contig_name();
}

inline GenomicRegion::SizeType get_begin(const GenomicRegion& a_region) noexcept
{
    return a_region.get_begin();
}

inline GenomicRegion::SizeType get_end(const GenomicRegion& a_region) noexcept
{
    return a_region.get_end();
}

inline std::string to_string(const GenomicRegion& a_region)
{
    return a_region.get_contig_name() + ':' + std::to_string(a_region.get_begin()) + '-'
    + std::to_string(a_region.get_end());
}

inline bool empty(const GenomicRegion& a_region) noexcept
{
    return empty(a_region.get_contig_region());
}

inline GenomicRegion::SizeType size(const GenomicRegion& a_region) noexcept
{
    return size(a_region.get_contig_region());
}

inline bool is_same_contig(const GenomicRegion& lhs, const GenomicRegion& rhs) noexcept
{
    return lhs.get_contig_name() == rhs.get_contig_name();
}

inline bool begins_equal(const GenomicRegion& lhs, const GenomicRegion& rhs)
{
    if (is_same_contig(lhs, rhs)) return begins_equal(lhs.get_contig_region(), rhs.get_contig_region());
    throw RegionError(to_string(lhs), to_string(rhs));
}

inline bool ends_equal(const GenomicRegion& lhs, const GenomicRegion& rhs)
{
    if (is_same_contig(lhs, rhs)) return ends_equal(lhs.get_contig_region(), rhs.get_contig_region());
    throw RegionError(to_string(lhs), to_string(rhs));
}

inline bool begins_before(const GenomicRegion& lhs, const GenomicRegion& rhs)
{
    if (is_same_contig(lhs, rhs)) return begins_before(lhs.get_contig_region(), rhs.get_contig_region());
    throw RegionError(to_string(lhs), to_string(rhs));
}

inline bool ends_before(const GenomicRegion& lhs, const GenomicRegion& rhs)
{
    if (is_same_contig(lhs, rhs)) return ends_before(lhs.get_contig_region(), rhs.get_contig_region());
    throw RegionError(to_string(lhs), to_string(rhs));
}

inline bool operator==(const GenomicRegion& lhs, const GenomicRegion& rhs) noexcept
{
    return is_same_contig(lhs, rhs) && lhs.get_contig_region() == rhs.get_contig_region();
}

inline bool operator<(const GenomicRegion& lhs, const GenomicRegion& rhs)
{
    if (is_same_contig(lhs, rhs)) return lhs.get_contig_region() < rhs.get_contig_region();
    throw RegionError(to_string(lhs), to_string(rhs));
}

inline bool is_before(const GenomicRegion& lhs, const GenomicRegion& rhs) noexcept
{
    if (is_same_contig(lhs, rhs)) return is_before(lhs.get_contig_region(), rhs.get_contig_region());
    throw RegionError(to_string(lhs), to_string(rhs));
}

inline bool is_after(const GenomicRegion& lhs, const GenomicRegion& rhs) noexcept
{
    if (is_same_contig(lhs, rhs)) return is_after(lhs.get_contig_region(), rhs.get_contig_region());
    throw RegionError(to_string(lhs), to_string(rhs));
}

inline bool are_adjacent(const GenomicRegion& lhs, const GenomicRegion& rhs) noexcept
{
    return is_same_contig(lhs, rhs) ? are_adjacent(lhs.get_contig_region(), rhs.get_contig_region()) : false;
}

inline GenomicRegion::DifferenceType overlap_size(const GenomicRegion& lhs, const GenomicRegion& rhs) noexcept
{
    return is_same_contig(lhs, rhs) ? overlap_size(lhs.get_contig_region(), rhs.get_contig_region()) : 0;
}

inline bool overlaps(const GenomicRegion& lhs, const GenomicRegion& rhs) noexcept
{
    return is_same_contig(lhs, rhs) ? overlaps(lhs.get_contig_region(), rhs.get_contig_region()) : false;
}

inline bool contains(const GenomicRegion& lhs, const GenomicRegion& rhs) noexcept
{
    return is_same_contig(lhs, rhs) && contains(lhs.get_contig_region(), rhs.get_contig_region());
}

inline GenomicRegion::DifferenceType inner_distance(const GenomicRegion& lhs, const GenomicRegion& rhs)
{
    if (is_same_contig(lhs, rhs)) {
        return inner_distance(lhs.get_contig_region(), rhs.get_contig_region());
    }
    throw RegionError(to_string(lhs), to_string(rhs));
}

inline GenomicRegion::DifferenceType outer_distance(const GenomicRegion& lhs, const GenomicRegion& rhs)
{
    if (is_same_contig(lhs, rhs)) {
        return outer_distance(lhs.get_contig_region(), rhs.get_contig_region());
    }
    throw RegionError(to_string(lhs), to_string(rhs));
}

inline GenomicRegion shift(const GenomicRegion& a_region, GenomicRegion::DifferenceType n)
{
    return GenomicRegion {
        a_region.get_contig_name(),
        shift(a_region.get_contig_region(), n)
    };
}

inline GenomicRegion next_position(const GenomicRegion& a_region)
{
    return GenomicRegion {a_region.get_contig_name(), next_position(a_region.get_contig_region())};
}

inline GenomicRegion compress_left(const GenomicRegion& a_region, GenomicRegion::DifferenceType n)
{
    return GenomicRegion {
        a_region.get_contig_name(),
        compress_left(a_region.get_contig_region(), n)
    };
}

inline GenomicRegion compress_right(const GenomicRegion& a_region, GenomicRegion::DifferenceType n)
{
    return GenomicRegion {
        a_region.get_contig_name(),
        compress_right(a_region.get_contig_region(), n)
    };
}

inline GenomicRegion get_encompassing(const GenomicRegion& lhs, const GenomicRegion& rhs)
{
    if (is_same_contig(lhs, rhs)) {
        return GenomicRegion {lhs.get_contig_name(), get_encompassing(lhs.get_contig_region(),
                                                                      rhs.get_contig_region())};
    }
    throw RegionError(to_string(lhs), to_string(rhs));
}

inline GenomicRegion get_intervening(const GenomicRegion& lhs, const GenomicRegion& rhs)
{
    return GenomicRegion {lhs.get_contig_name(), get_intervening(lhs.get_contig_region(),
                                                                 rhs.get_contig_region())};
}

inline GenomicRegion get_overlapped(const GenomicRegion& lhs, const GenomicRegion& rhs) noexcept
{
    return GenomicRegion {lhs.get_contig_name(), get_overlapped(lhs.get_contig_region(),
                                                                rhs.get_contig_region())};
}

inline GenomicRegion get_left_overhang(const GenomicRegion& lhs, const GenomicRegion& rhs)
{
    if (is_same_contig(lhs, rhs)) {
        return GenomicRegion {lhs.get_contig_name(), get_left_overhang(lhs.get_contig_region(),
                                                                       rhs.get_contig_region())};
    }
    throw RegionError(to_string(lhs), to_string(rhs));
}

inline GenomicRegion get_right_overhang(const GenomicRegion& lhs, const GenomicRegion& rhs)
{
    if (is_same_contig(lhs, rhs)) {
        return GenomicRegion {lhs.get_contig_name(), get_right_overhang(lhs.get_contig_region(),
                                                                        rhs.get_contig_region())};
    }
    throw RegionError(to_string(lhs), to_string(rhs));
}

inline GenomicRegion get_closed(const GenomicRegion& lhs, const GenomicRegion& rhs)
{
    if (is_same_contig(lhs, rhs)) {
        return GenomicRegion {lhs.get_contig_name(), get_closed(lhs.get_contig_region(),
                                                                rhs.get_contig_region())};
    }
    throw RegionError(to_string(lhs), to_string(rhs));
}

namespace std {
    template <> struct hash<GenomicRegion>
    {
        size_t operator()(const GenomicRegion& r) const
        {
            size_t seed {};
            boost::hash_combine(seed, hash<GenomicRegion::StringType>()(r.get_contig_name()));
            boost::hash_combine(seed, hash<ContigRegion>()(r.get_contig_region()));
            return seed;
        }
    };
}

inline std::ostream& operator<<(std::ostream& os, const GenomicRegion& a_region)
{
    os << to_string(a_region);
    return os;
}


#endif /* defined(__Octopus__genomic_region__) */
