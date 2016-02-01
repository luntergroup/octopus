//
//  mappable.hpp
//  Octopus
//
//  Created by Daniel Cooke on 10/04/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_mappable_hpp
#define Octopus_mappable_hpp

#include <type_traits>

#include "contig_region.hpp"
#include "genomic_region.hpp"

/**
 Template black magic. Anything that inherits from Mappable and has a get_region() member method
 can use all of these methods.
*/

template <typename T>
class Mappable {};

template <typename T>
constexpr bool is_contig_region = std::is_same<std::decay_t<T>, ContigRegion>::value;

template <typename T>
constexpr bool is_genomic_region = std::is_same<std::decay_t<T>, GenomicRegion>::value;

template <typename T>
constexpr bool is_region = is_contig_region<T> || is_genomic_region<T>;

template <typename T>
constexpr bool is_mappable = std::is_base_of<Mappable<T>, T>::value;

template <typename T>
constexpr bool is_region_or_mappable = is_region<T> || is_mappable<T>;

template <typename T, typename R = void>
using EnableIfMappable = std::enable_if_t<is_mappable<T>, R>;

template <typename T, typename R = void>
using EnableIfRegionOrMappable = std::enable_if_t<is_region_or_mappable<T>, R>;

inline decltype(auto) get_region(const ContigRegion& region)
{
    return region;
}

inline decltype(auto) get_region(const GenomicRegion& region)
{
    return region;
}

template <typename T>
inline decltype(auto) get_region(const Mappable<T>& m)
{
    return static_cast<const T&>(m).get_region();
}

template <typename T>
using RegionType = std::decay_t<decltype(get_region(std::declval<std::decay_t<T>>()))>;

inline decltype(auto) get_contig_region(const ContigRegion& region) noexcept
{
    return region;
}

inline decltype(auto) get_contig_region(const GenomicRegion& region) noexcept
{
    return region.get_contig_region();
}

template <typename T>
inline decltype(auto) get_contig_region(const T& mappable) noexcept
{
    return get_contig_region(static_cast<const T&>(mappable).get_region());
}

inline auto is_same_region(const ContigRegion& lhs, const ContigRegion& rhs)
{
    return lhs == rhs;
}

inline auto is_same_region(const GenomicRegion& lhs, const GenomicRegion& rhs)
{
    return lhs == rhs;
}

template <typename T>
inline auto is_same_region(const ContigRegion& lhs, const Mappable<T>& rhs)
{
    return is_same_region(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T>
inline auto is_same_region(const GenomicRegion& lhs, const Mappable<T>& rhs)
{
    return is_same_region(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T>
inline auto is_same_region(const Mappable<T>& lhs, const ContigRegion& rhs)
{
    return is_same_region(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
inline auto is_same_region(const Mappable<T>& lhs, const GenomicRegion& rhs)
{
    return is_same_region(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T1, typename T2>
inline auto is_same_region(const Mappable<T1>& lhs, const Mappable<T2>& rhs)
{
    return is_same_region(static_cast<const T1&>(lhs).get_region(),
                          static_cast<const T2&>(rhs).get_region());
}

template <typename T>
inline auto empty(const Mappable<T>& m)
{
    return empty(static_cast<const T&>(m).get_region());
}

template <typename T>
inline auto size(const Mappable<T>& m)
{
    return size(static_cast<const T&>(m).get_region());
}

template <typename T>
inline decltype(auto) get_contig_name(const Mappable<T>& m)
{
    return get_contig_name(static_cast<const T&>(m).get_region());
}

template <typename T>
inline auto is_same_contig(const Mappable<T>& lhs, const ContigRegion& rhs)
{
    return is_same_contig(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
inline auto is_same_contig(const Mappable<T>& lhs, const GenomicRegion& rhs)
{
    return is_same_contig(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T1, typename T2>
inline auto is_same_contig(const Mappable<T1>& lhs, const Mappable<T2>& rhs)
{
    return is_same_contig(static_cast<const T1&>(lhs).get_region(),
                          static_cast<const T2&>(rhs).get_region());
}

template <typename T>
inline auto get_begin(const Mappable<T>& m)
{
    return get_begin(static_cast<const T&>(m).get_region());
}

template <typename T>
inline auto get_end(const Mappable<T>& m)
{
    return get_end(static_cast<const T&>(m).get_region());
}

template <typename T>
inline auto operator==(const Mappable<T>& lhs, const ContigRegion& rhs)
{
    return static_cast<const T&>(lhs).get_region() == rhs;
}

template <typename T>
inline auto operator==(const Mappable<T>& lhs, const GenomicRegion& rhs)
{
    return static_cast<const T&>(lhs).get_region() == rhs;
}

template <typename T>
inline auto operator==(const ContigRegion& lhs, const Mappable<T>& rhs)
{
    return lhs == static_cast<const T&>(rhs).get_region();
}

template <typename T>
inline auto operator==(const GenomicRegion& lhs, const Mappable<T>& rhs)
{
    return lhs == static_cast<const T&>(rhs).get_region();
}

template <typename T1, typename T2>
inline auto operator==(const Mappable<T1>& lhs, const Mappable<T2>& rhs)
{
    return static_cast<const T1&>(lhs).get_region() == static_cast<const T2&>(rhs).get_region();
}

template <typename T>
inline auto operator<(const Mappable<T>& lhs, const ContigRegion& rhs)
{
    return static_cast<const T&>(lhs).get_region() < rhs;
}

template <typename T>
inline auto operator<(const Mappable<T>& lhs, const GenomicRegion& rhs)
{
    return static_cast<const T&>(lhs).get_region() < rhs;
}

template <typename T>
inline auto operator<(const ContigRegion& lhs, const Mappable<T>& rhs)
{
    return lhs < static_cast<const T&>(rhs).get_region();
}

template <typename T>
inline auto operator<(const GenomicRegion& lhs, const Mappable<T>& rhs)
{
    return lhs < static_cast<const T&>(rhs).get_region();
}

template <typename T1, typename T2>
inline auto operator<(const Mappable<T1>& lhs, const Mappable<T2>& rhs)
{
    return static_cast<const T1&>(lhs).get_region() < static_cast<const T2&>(rhs).get_region();
}

template <typename T>
inline auto begins_equal(const ContigRegion& lhs, const Mappable<T>& rhs)
{
    return begins_equal(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T>
inline auto begins_equal(const GenomicRegion& lhs, const Mappable<T>& rhs)
{
    return begins_equal(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T>
inline auto begins_equal(const Mappable<T>& lhs, const ContigRegion& rhs)
{
    return begins_equal(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
inline auto begins_equal(const Mappable<T>& lhs, const GenomicRegion& rhs)
{
    return begins_equal(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T1, typename T2>
inline auto begins_equal(const Mappable<T1>& lhs, const Mappable<T2>& rhs)
{
    return begins_equal(static_cast<const T1&>(lhs).get_region(),
                        static_cast<const T2&>(rhs).get_region());
}

template <typename T>
inline auto ends_equal(const ContigRegion& lhs, const Mappable<T>& rhs)
{
    return ends_equal(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T>
inline auto ends_equal(const GenomicRegion& lhs, const Mappable<T>& rhs)
{
    return ends_equal(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T>
inline auto ends_equal(const Mappable<T>& lhs, const ContigRegion& rhs)
{
    return ends_equal(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
inline auto ends_equal(const Mappable<T>& lhs, const GenomicRegion& rhs)
{
    return ends_equal(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T1, typename T2>
inline auto ends_equal(const Mappable<T1>& lhs, const Mappable<T2>& rhs)
{
    return ends_equal(static_cast<const T1&>(lhs).get_region(),
                      static_cast<const T2&>(rhs).get_region());
}

template <typename T>
inline auto begins_before(const Mappable<T>& lhs, const ContigRegion& rhs)
{
    return begins_before(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
inline auto begins_before(const Mappable<T>& lhs, const GenomicRegion& rhs)
{
    return begins_before(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
inline auto begins_before(const ContigRegion& lhs, const Mappable<T>& rhs)
{
    return begins_before(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T>
inline auto begins_before(const GenomicRegion& lhs, const Mappable<T>& rhs)
{
    return begins_before(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T1, typename T2>
inline auto begins_before(const Mappable<T1>& lhs, const Mappable<T2>& rhs)
{
    return begins_before(static_cast<const T1&>(lhs).get_region(),
                         static_cast<const T2&>(rhs).get_region());
}

template <typename T>
inline auto ends_before(const Mappable<T>& lhs, const ContigRegion& rhs)
{
    return ends_before(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
inline auto ends_before(const Mappable<T>& lhs, const GenomicRegion& rhs)
{
    return ends_before(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
inline auto ends_before(const ContigRegion& lhs, const Mappable<T>& rhs)
{
    return ends_before(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T>
inline auto ends_before(const GenomicRegion& lhs, const Mappable<T>& rhs)
{
    return ends_before(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T1, typename T2>
inline auto ends_before(const Mappable<T1>& lhs, const Mappable<T2>& rhs)
{
    return ends_before(static_cast<const T1&>(lhs).get_region(),
                       static_cast<const T2&>(rhs).get_region());
}

template <typename T>
inline auto is_before(const Mappable<T>& lhs, const ContigRegion& rhs) noexcept
{
    return is_before(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
inline auto is_before(const Mappable<T>& lhs, const GenomicRegion& rhs) noexcept
{
    return is_before(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
inline auto is_before(const ContigRegion& lhs, const Mappable<T>& rhs) noexcept
{
    return is_before(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T>
inline auto is_before(const GenomicRegion& lhs, const Mappable<T>& rhs) noexcept
{
    return is_before(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T1, typename T2>
inline auto is_before(const Mappable<T1>& lhs, const Mappable<T2>& rhs) noexcept
{
    return is_before(static_cast<const T1&>(lhs).get_region(),
                     static_cast<const T2&>(rhs).get_region());
}

template <typename T>
inline auto is_after(const Mappable<T>& lhs, const ContigRegion& rhs) noexcept
{
    return is_after(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
inline auto is_after(const Mappable<T>& lhs, const GenomicRegion& rhs) noexcept
{
    return is_after(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
inline auto is_after(const ContigRegion& lhs, const Mappable<T>& rhs) noexcept
{
    return is_after(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T>
inline auto is_after(const GenomicRegion& lhs, const Mappable<T>& rhs) noexcept
{
    return is_after(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T1, typename T2>
inline auto is_after(const Mappable<T1>& lhs, const Mappable<T2>& rhs) noexcept
{
    return is_after(static_cast<const T1&>(lhs).get_region(),
                    static_cast<const T2&>(rhs).get_region());
}

template <typename T>
inline auto are_adjacent(const Mappable<T>& lhs, const ContigRegion& rhs) noexcept
{
    return are_adjacent(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
inline auto are_adjacent(const Mappable<T>& lhs, const GenomicRegion& rhs) noexcept
{
    return are_adjacent(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
inline auto are_adjacent(const ContigRegion& lhs, const Mappable<T>& rhs) noexcept
{
    return are_adjacent(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T>
inline auto are_adjacent(const GenomicRegion& lhs, const Mappable<T>& rhs) noexcept
{
    return are_adjacent(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T1, typename T2>
inline auto are_adjacent(const Mappable<T1>& lhs, const Mappable<T2>& rhs) noexcept
{
    return are_adjacent(static_cast<const T1&>(lhs).get_region(),
                        static_cast<const T2&>(rhs).get_region());
}

template <typename T>
inline auto overlap_size(const Mappable<T>& lhs, const GenomicRegion& rhs) noexcept
{
    return overlap_size(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
inline auto overlap_size(const GenomicRegion& lhs, const Mappable<T>& rhs) noexcept
{
    return overlap_size(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T1, typename T2>
inline auto overlap_size(const Mappable<T1>& lhs, const Mappable<T2>& rhs) noexcept
{
    return overlap_size(static_cast<const T1&>(lhs).get_region(),
                        static_cast<const T2&>(rhs).get_region());
}

template <typename T>
inline auto overlaps(const Mappable<T>& lhs, const ContigRegion& rhs) noexcept
{
    return overlaps(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
inline auto overlaps(const Mappable<T>& lhs, const GenomicRegion& rhs) noexcept
{
    return overlaps(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
inline auto overlaps(const ContigRegion& lhs, const Mappable<T>& rhs) noexcept
{
    return overlaps(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T>
inline auto overlaps(const GenomicRegion& lhs, const Mappable<T>& rhs) noexcept
{
    return overlaps(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T1, typename T2>
inline auto overlaps(const Mappable<T1>& lhs, const Mappable<T2>& rhs) noexcept
{
    return overlaps(static_cast<const T1&>(lhs).get_region(),
                    static_cast<const T2&>(rhs).get_region());
}

template <typename T>
inline auto contains(const Mappable<T>& lhs, const ContigRegion& rhs) noexcept
{
    return contains(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
inline auto contains(const Mappable<T>& lhs, const GenomicRegion& rhs) noexcept
{
    return contains(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
inline auto contains(const ContigRegion& lhs, const Mappable<T>& rhs) noexcept
{
    return contains(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T>
inline auto contains(const GenomicRegion& lhs, const Mappable<T>& rhs) noexcept
{
    return contains(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T1, typename T2>
inline auto contains(const Mappable<T1>& lhs, const Mappable<T2>& rhs) noexcept
{
    return contains(static_cast<const T1&>(lhs).get_region(),
                    static_cast<const T2&>(rhs).get_region());
}

template <typename T>
inline auto inner_distance(const Mappable<T>& lhs, const ContigRegion& rhs)
{
    return inner_distance(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
inline auto inner_distance(const Mappable<T>& lhs, const GenomicRegion& rhs)
{
    return inner_distance(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
inline auto inner_distance(const ContigRegion& lhs, const Mappable<T>& rhs)
{
    return inner_distance(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T>
inline auto inner_distance(const GenomicRegion& lhs, const Mappable<T>& rhs)
{
    return inner_distance(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T1, typename T2>
inline auto inner_distance(const Mappable<T1>& lhs, const Mappable<T2>& rhs)
{
    return inner_distance(static_cast<const T1&>(lhs).get_region(),
                          static_cast<const T2&>(rhs).get_region());
}

template <typename T>
inline auto outer_distance(const Mappable<T>& lhs, const ContigRegion& rhs)
{
    return outer_distance(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
inline auto outer_distance(const Mappable<T>& lhs, const GenomicRegion& rhs)
{
    return outer_distance(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
inline auto outer_distance(const ContigRegion& lhs, const Mappable<T>& rhs)
{
    return outer_distance(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T>
inline auto outer_distance(const GenomicRegion& lhs, const Mappable<T>& rhs)
{
    return outer_distance(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T1, typename T2>
inline auto outer_distance(const Mappable<T1>& lhs, const Mappable<T2>& rhs)
{
    return outer_distance(static_cast<const T1&>(lhs).get_region(),
                          static_cast<const T2&>(rhs).get_region());
}

template <typename T>
inline auto next_position(const Mappable<T>& mappable)
{
    return next_position(static_cast<const T&>(mappable).get_region());
}

template <typename T>
inline auto compress_lhs(const Mappable<T>& mappable, const GenomicRegion::DifferenceType n)
{
    return compress_lhs(static_cast<const T&>(mappable).get_region(), n);
}

template <typename T>
inline auto compress_rhs(const Mappable<T>& mappable, const GenomicRegion::DifferenceType n)
{
    return compress_rhs(static_cast<const T&>(mappable).get_region(), n);
}

template <typename T>
inline auto get_encompassing(const Mappable<T>& lhs, const ContigRegion& rhs)
{
    return get_encompassing(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
inline auto get_encompassing(const Mappable<T>& lhs, const GenomicRegion& rhs)
{
    return get_encompassing(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
inline auto get_encompassing(const ContigRegion& lhs, const Mappable<T>& rhs)
{
    return get_encompassing(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T>
inline auto get_encompassing(const GenomicRegion& lhs, const Mappable<T>& rhs)
{
    return get_encompassing(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T1, typename T2>
inline auto get_encompassing(const Mappable<T1>& lhs, const Mappable<T2>& rhs)
{
    return get_encompassing(static_cast<const T1&>(lhs).get_region(),
                            static_cast<const T2&>(rhs).get_region());
}

template <typename T>
inline auto get_intervening(const Mappable<T>& lhs, const ContigRegion& rhs)
{
    return get_intervening(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
inline auto get_intervening(const Mappable<T>& lhs, const GenomicRegion& rhs)
{
    return get_intervening(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
inline auto get_intervening(const ContigRegion& lhs, const Mappable<T>& rhs)
{
    return get_intervening(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T>
inline auto get_intervening(const GenomicRegion& lhs, const Mappable<T>& rhs)
{
    return get_intervening(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T1, typename T2>
inline auto get_intervening(const Mappable<T1>& lhs, const Mappable<T2>& rhs)
{
    return get_intervening(static_cast<const T1&>(lhs).get_region(),
                           static_cast<const T2&>(rhs).get_region());
}


template <typename T>
inline auto get_overlapped(const Mappable<T>& lhs, const ContigRegion& rhs)
{
    return get_overlapped(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
inline auto get_overlapped(const Mappable<T>& lhs, const GenomicRegion& rhs)
{
    return get_overlapped(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
inline auto get_overlapped(const ContigRegion& lhs, const Mappable<T>& rhs)
{
    return get_overlapped(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T>
inline auto get_overlapped(const GenomicRegion& lhs, const Mappable<T>& rhs)
{
    return get_overlapped(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T1, typename T2>
inline auto get_overlapped(const Mappable<T1>& lhs, const Mappable<T2>& rhs)
{
    return get_overlapped(static_cast<const T1&>(lhs).get_region(),
                          static_cast<const T2&>(rhs).get_region());
}

template <typename T>
inline auto left_overhang_size(const Mappable<T>& lhs, const ContigRegion& rhs)
{
    return left_overhang_size(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
inline auto left_overhang_size(const Mappable<T>& lhs, const GenomicRegion& rhs)
{
    return left_overhang_size(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
inline auto left_overhang_size(const ContigRegion& lhs, const Mappable<T>& rhs)
{
    return left_overhang_size(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T>
inline auto left_overhang_size(const GenomicRegion& lhs, const Mappable<T>& rhs)
{
    return left_overhang_size(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T1, typename T2>
inline auto left_overhang_size(const Mappable<T1>& lhs, const Mappable<T2>& rhs)
{
    return left_overhang_size(static_cast<const T1&>(lhs).get_region(),
                              static_cast<const T2&>(rhs).get_region());
}

template <typename T>
inline auto right_overhang_size(const Mappable<T>& lhs, const ContigRegion& rhs)
{
    return right_overhang_size(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
inline auto right_overhang_size(const Mappable<T>& lhs, const GenomicRegion& rhs)
{
    return right_overhang_size(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
inline auto right_overhang_size(const ContigRegion& lhs, const Mappable<T>& rhs)
{
    return right_overhang_size(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T>
inline auto right_overhang_size(const GenomicRegion& lhs, const Mappable<T>& rhs)
{
    return right_overhang_size(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T1, typename T2>
inline auto right_overhang_size(const Mappable<T1>& lhs, const Mappable<T2>& rhs)
{
    return right_overhang_size(static_cast<const T1&>(lhs).get_region(),
                               static_cast<const T2&>(rhs).get_region());
}

template <typename T>
inline auto get_left_overhang(const Mappable<T>& lhs, const ContigRegion& rhs)
{
    return get_left_overhang(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
inline auto get_left_overhang(const Mappable<T>& lhs, const GenomicRegion& rhs)
{
    return get_left_overhang(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
inline auto get_left_overhang(const ContigRegion& lhs, const Mappable<T>& rhs)
{
    return get_left_overhang(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T>
inline auto get_left_overhang(const GenomicRegion& lhs, const Mappable<T>& rhs)
{
    return get_left_overhang(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T1, typename T2>
inline auto get_left_overhang(const Mappable<T1>& lhs, const Mappable<T2>& rhs)
{
    return get_left_overhang(static_cast<const T1&>(lhs).get_region(),
                             static_cast<const T2&>(rhs).get_region());
}

template <typename T>
inline auto get_right_overhang(const Mappable<T>& lhs, const ContigRegion& rhs)
{
    return get_right_overhang(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
inline auto get_right_overhang(const Mappable<T>& lhs, const GenomicRegion& rhs)
{
    return get_right_overhang(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
inline auto get_right_overhang(const ContigRegion& lhs, const Mappable<T>& rhs)
{
    return get_right_overhang(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T>
inline auto get_right_overhang(const GenomicRegion& lhs, const Mappable<T>& rhs)
{
    return get_right_overhang(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T1, typename T2>
inline auto get_right_overhang(const Mappable<T1>& lhs, const Mappable<T2>& rhs)
{
    return get_right_overhang(static_cast<const T1&>(lhs).get_region(),
                              static_cast<const T2&>(rhs).get_region());
}

template <typename T>
inline auto get_closed(const Mappable<T>& lhs, const ContigRegion& rhs)
{
    return get_closed(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
inline auto get_closed(const Mappable<T>& lhs, const GenomicRegion& rhs)
{
    return get_closed(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
inline auto get_closed(const ContigRegion& lhs, const Mappable<T>& rhs)
{
    return get_closed(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T>
inline auto get_closed(const GenomicRegion& lhs, const Mappable<T>& rhs)
{
    return get_closed(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T1, typename T2>
inline auto get_closed(const Mappable<T1>& lhs, const Mappable<T2>& rhs)
{
    return get_closed(static_cast<const T1&>(lhs).get_region(),
                      static_cast<const T2&>(rhs).get_region());
}

template <typename T>
inline auto get_head(const Mappable<T>& m, const GenomicRegion::SizeType n = 0)
{
    return get_head(static_cast<const T&>(m).get_region(), n);
}

template <typename T>
inline auto get_tail(const Mappable<T>& m, const GenomicRegion::SizeType n = 0)
{
    return get_tail(static_cast<const T&>(m).get_region(), n);
}

#endif
