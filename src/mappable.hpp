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
decltype(auto) get_region(const Mappable<T>& m)
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
decltype(auto) get_contig_region(const T& mappable) noexcept
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
auto is_same_region(const ContigRegion& lhs, const Mappable<T>& rhs)
{
    return is_same_region(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T>
auto is_same_region(const GenomicRegion& lhs, const Mappable<T>& rhs)
{
    return is_same_region(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T>
auto is_same_region(const Mappable<T>& lhs, const ContigRegion& rhs)
{
    return is_same_region(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
auto is_same_region(const Mappable<T>& lhs, const GenomicRegion& rhs)
{
    return is_same_region(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T1, typename T2>
auto is_same_region(const Mappable<T1>& lhs, const Mappable<T2>& rhs)
{
    return is_same_region(static_cast<const T1&>(lhs).get_region(),
                          static_cast<const T2&>(rhs).get_region());
}

template <typename T>
auto empty(const Mappable<T>& m)
{
    return empty(static_cast<const T&>(m).get_region());
}

template <typename T>
auto size(const Mappable<T>& m)
{
    return size(static_cast<const T&>(m).get_region());
}

template <typename T>
decltype(auto) get_contig_name(const Mappable<T>& m)
{
    return get_contig_name(static_cast<const T&>(m).get_region());
}

template <typename T>
auto is_same_contig(const Mappable<T>& lhs, const ContigRegion& rhs)
{
    return is_same_contig(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
auto is_same_contig(const Mappable<T>& lhs, const GenomicRegion& rhs)
{
    return is_same_contig(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T1, typename T2>
auto is_same_contig(const Mappable<T1>& lhs, const Mappable<T2>& rhs)
{
    return is_same_contig(static_cast<const T1&>(lhs).get_region(),
                          static_cast<const T2&>(rhs).get_region());
}

template <typename T>
auto get_begin(const Mappable<T>& m)
{
    return get_begin(static_cast<const T&>(m).get_region());
}

template <typename T>
auto get_end(const Mappable<T>& m)
{
    return get_end(static_cast<const T&>(m).get_region());
}

template <typename T>
auto operator==(const Mappable<T>& lhs, const ContigRegion& rhs)
{
    return static_cast<const T&>(lhs).get_region() == rhs;
}

template <typename T>
auto operator==(const Mappable<T>& lhs, const GenomicRegion& rhs)
{
    return static_cast<const T&>(lhs).get_region() == rhs;
}

template <typename T>
auto operator==(const ContigRegion& lhs, const Mappable<T>& rhs)
{
    return lhs == static_cast<const T&>(rhs).get_region();
}

template <typename T>
auto operator==(const GenomicRegion& lhs, const Mappable<T>& rhs)
{
    return lhs == static_cast<const T&>(rhs).get_region();
}

template <typename T1, typename T2>
auto operator==(const Mappable<T1>& lhs, const Mappable<T2>& rhs)
{
    return static_cast<const T1&>(lhs).get_region() == static_cast<const T2&>(rhs).get_region();
}

template <typename T>
auto operator<(const Mappable<T>& lhs, const ContigRegion& rhs)
{
    return static_cast<const T&>(lhs).get_region() < rhs;
}

template <typename T>
auto operator<(const Mappable<T>& lhs, const GenomicRegion& rhs)
{
    return static_cast<const T&>(lhs).get_region() < rhs;
}

template <typename T>
auto operator<(const ContigRegion& lhs, const Mappable<T>& rhs)
{
    return lhs < static_cast<const T&>(rhs).get_region();
}

template <typename T>
auto operator<(const GenomicRegion& lhs, const Mappable<T>& rhs)
{
    return lhs < static_cast<const T&>(rhs).get_region();
}

template <typename T1, typename T2>
auto operator<(const Mappable<T1>& lhs, const Mappable<T2>& rhs)
{
    return static_cast<const T1&>(lhs).get_region() < static_cast<const T2&>(rhs).get_region();
}

template <typename T>
auto begins_equal(const ContigRegion& lhs, const Mappable<T>& rhs)
{
    return begins_equal(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T>
auto begins_equal(const GenomicRegion& lhs, const Mappable<T>& rhs)
{
    return begins_equal(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T>
auto begins_equal(const Mappable<T>& lhs, const ContigRegion& rhs)
{
    return begins_equal(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
auto begins_equal(const Mappable<T>& lhs, const GenomicRegion& rhs)
{
    return begins_equal(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T1, typename T2>
auto begins_equal(const Mappable<T1>& lhs, const Mappable<T2>& rhs)
{
    return begins_equal(static_cast<const T1&>(lhs).get_region(),
                        static_cast<const T2&>(rhs).get_region());
}

template <typename T>
auto ends_equal(const ContigRegion& lhs, const Mappable<T>& rhs)
{
    return ends_equal(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T>
auto ends_equal(const GenomicRegion& lhs, const Mappable<T>& rhs)
{
    return ends_equal(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T>
auto ends_equal(const Mappable<T>& lhs, const ContigRegion& rhs)
{
    return ends_equal(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
auto ends_equal(const Mappable<T>& lhs, const GenomicRegion& rhs)
{
    return ends_equal(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T1, typename T2>
auto ends_equal(const Mappable<T1>& lhs, const Mappable<T2>& rhs)
{
    return ends_equal(static_cast<const T1&>(lhs).get_region(),
                      static_cast<const T2&>(rhs).get_region());
}

template <typename T>
auto begins_before(const Mappable<T>& lhs, const ContigRegion& rhs)
{
    return begins_before(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
auto begins_before(const Mappable<T>& lhs, const GenomicRegion& rhs)
{
    return begins_before(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
auto begins_before(const ContigRegion& lhs, const Mappable<T>& rhs)
{
    return begins_before(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T>
auto begins_before(const GenomicRegion& lhs, const Mappable<T>& rhs)
{
    return begins_before(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T1, typename T2>
auto begins_before(const Mappable<T1>& lhs, const Mappable<T2>& rhs)
{
    return begins_before(static_cast<const T1&>(lhs).get_region(),
                         static_cast<const T2&>(rhs).get_region());
}

template <typename T>
auto ends_before(const Mappable<T>& lhs, const ContigRegion& rhs)
{
    return ends_before(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
auto ends_before(const Mappable<T>& lhs, const GenomicRegion& rhs)
{
    return ends_before(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
auto ends_before(const ContigRegion& lhs, const Mappable<T>& rhs)
{
    return ends_before(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T>
auto ends_before(const GenomicRegion& lhs, const Mappable<T>& rhs)
{
    return ends_before(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T1, typename T2>
auto ends_before(const Mappable<T1>& lhs, const Mappable<T2>& rhs)
{
    return ends_before(static_cast<const T1&>(lhs).get_region(),
                       static_cast<const T2&>(rhs).get_region());
}

template <typename T>
auto is_before(const Mappable<T>& lhs, const ContigRegion& rhs) noexcept
{
    return is_before(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
auto is_before(const Mappable<T>& lhs, const GenomicRegion& rhs) noexcept
{
    return is_before(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
auto is_before(const ContigRegion& lhs, const Mappable<T>& rhs) noexcept
{
    return is_before(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T>
auto is_before(const GenomicRegion& lhs, const Mappable<T>& rhs) noexcept
{
    return is_before(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T1, typename T2>
auto is_before(const Mappable<T1>& lhs, const Mappable<T2>& rhs) noexcept
{
    return is_before(static_cast<const T1&>(lhs).get_region(),
                     static_cast<const T2&>(rhs).get_region());
}

template <typename T>
auto is_after(const Mappable<T>& lhs, const ContigRegion& rhs) noexcept
{
    return is_after(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
auto is_after(const Mappable<T>& lhs, const GenomicRegion& rhs) noexcept
{
    return is_after(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
auto is_after(const ContigRegion& lhs, const Mappable<T>& rhs) noexcept
{
    return is_after(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T>
auto is_after(const GenomicRegion& lhs, const Mappable<T>& rhs) noexcept
{
    return is_after(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T1, typename T2>
auto is_after(const Mappable<T1>& lhs, const Mappable<T2>& rhs) noexcept
{
    return is_after(static_cast<const T1&>(lhs).get_region(),
                    static_cast<const T2&>(rhs).get_region());
}

template <typename T>
auto are_adjacent(const Mappable<T>& lhs, const ContigRegion& rhs) noexcept
{
    return are_adjacent(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
auto are_adjacent(const Mappable<T>& lhs, const GenomicRegion& rhs) noexcept
{
    return are_adjacent(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
auto are_adjacent(const ContigRegion& lhs, const Mappable<T>& rhs) noexcept
{
    return are_adjacent(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T>
auto are_adjacent(const GenomicRegion& lhs, const Mappable<T>& rhs) noexcept
{
    return are_adjacent(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T1, typename T2>
auto are_adjacent(const Mappable<T1>& lhs, const Mappable<T2>& rhs) noexcept
{
    return are_adjacent(static_cast<const T1&>(lhs).get_region(),
                        static_cast<const T2&>(rhs).get_region());
}

template <typename T>
auto overlap_size(const Mappable<T>& lhs, const GenomicRegion& rhs) noexcept
{
    return overlap_size(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
auto overlap_size(const GenomicRegion& lhs, const Mappable<T>& rhs) noexcept
{
    return overlap_size(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T1, typename T2>
auto overlap_size(const Mappable<T1>& lhs, const Mappable<T2>& rhs) noexcept
{
    return overlap_size(static_cast<const T1&>(lhs).get_region(),
                        static_cast<const T2&>(rhs).get_region());
}

template <typename T>
auto overlaps(const Mappable<T>& lhs, const ContigRegion& rhs) noexcept
{
    return overlaps(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
auto overlaps(const Mappable<T>& lhs, const GenomicRegion& rhs) noexcept
{
    return overlaps(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
auto overlaps(const ContigRegion& lhs, const Mappable<T>& rhs) noexcept
{
    return overlaps(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T>
auto overlaps(const GenomicRegion& lhs, const Mappable<T>& rhs) noexcept
{
    return overlaps(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T1, typename T2>
auto overlaps(const Mappable<T1>& lhs, const Mappable<T2>& rhs) noexcept
{
    return overlaps(static_cast<const T1&>(lhs).get_region(),
                    static_cast<const T2&>(rhs).get_region());
}

template <typename T>
auto contains(const Mappable<T>& lhs, const ContigRegion& rhs) noexcept
{
    return contains(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
auto contains(const Mappable<T>& lhs, const GenomicRegion& rhs) noexcept
{
    return contains(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
auto contains(const ContigRegion& lhs, const Mappable<T>& rhs) noexcept
{
    return contains(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T>
auto contains(const GenomicRegion& lhs, const Mappable<T>& rhs) noexcept
{
    return contains(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T1, typename T2>
auto contains(const Mappable<T1>& lhs, const Mappable<T2>& rhs) noexcept
{
    return contains(static_cast<const T1&>(lhs).get_region(),
                    static_cast<const T2&>(rhs).get_region());
}

template <typename T>
auto inner_distance(const Mappable<T>& lhs, const ContigRegion& rhs)
{
    return inner_distance(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
auto inner_distance(const Mappable<T>& lhs, const GenomicRegion& rhs)
{
    return inner_distance(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
auto inner_distance(const ContigRegion& lhs, const Mappable<T>& rhs)
{
    return inner_distance(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T>
auto inner_distance(const GenomicRegion& lhs, const Mappable<T>& rhs)
{
    return inner_distance(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T1, typename T2>
auto inner_distance(const Mappable<T1>& lhs, const Mappable<T2>& rhs)
{
    return inner_distance(static_cast<const T1&>(lhs).get_region(),
                          static_cast<const T2&>(rhs).get_region());
}

template <typename T>
auto outer_distance(const Mappable<T>& lhs, const ContigRegion& rhs)
{
    return outer_distance(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
auto outer_distance(const Mappable<T>& lhs, const GenomicRegion& rhs)
{
    return outer_distance(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
auto outer_distance(const ContigRegion& lhs, const Mappable<T>& rhs)
{
    return outer_distance(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T>
auto outer_distance(const GenomicRegion& lhs, const Mappable<T>& rhs)
{
    return outer_distance(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T1, typename T2>
auto outer_distance(const Mappable<T1>& lhs, const Mappable<T2>& rhs)
{
    return outer_distance(static_cast<const T1&>(lhs).get_region(),
                          static_cast<const T2&>(rhs).get_region());
}

template <typename T>
auto next_position(const Mappable<T>& mappable)
{
    return next_position(static_cast<const T&>(mappable).get_region());
}

template <typename T>
auto compress_lhs(const Mappable<T>& mappable, const GenomicRegion::DifferenceType n)
{
    return compress_lhs(static_cast<const T&>(mappable).get_region(), n);
}

template <typename T>
auto compress_rhs(const Mappable<T>& mappable, const GenomicRegion::DifferenceType n)
{
    return compress_rhs(static_cast<const T&>(mappable).get_region(), n);
}

template <typename T>
auto encompassing_region(const Mappable<T>& lhs, const ContigRegion& rhs)
{
    return encompassing_region(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
auto encompassing_region(const Mappable<T>& lhs, const GenomicRegion& rhs)
{
    return encompassing_region(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
auto encompassing_region(const ContigRegion& lhs, const Mappable<T>& rhs)
{
    return encompassing_region(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T>
auto encompassing_region(const GenomicRegion& lhs, const Mappable<T>& rhs)
{
    return encompassing_region(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T1, typename T2>
auto encompassing_region(const Mappable<T1>& lhs, const Mappable<T2>& rhs)
{
    return encompassing_region(static_cast<const T1&>(lhs).get_region(),
                               static_cast<const T2&>(rhs).get_region());
}

template <typename T>
auto intervening_region(const Mappable<T>& lhs, const ContigRegion& rhs)
{
    return intervening_region(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
auto intervening_region(const Mappable<T>& lhs, const GenomicRegion& rhs)
{
    return intervening_region(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
auto intervening_region(const ContigRegion& lhs, const Mappable<T>& rhs)
{
    return intervening_region(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T>
auto intervening_region(const GenomicRegion& lhs, const Mappable<T>& rhs)
{
    return intervening_region(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T1, typename T2>
auto intervening_region(const Mappable<T1>& lhs, const Mappable<T2>& rhs)
{
    return intervening_region(static_cast<const T1&>(lhs).get_region(),
                              static_cast<const T2&>(rhs).get_region());
}


template <typename T>
auto overlapped_region(const Mappable<T>& lhs, const ContigRegion& rhs)
{
    return overlapped_region(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
auto overlapped_region(const Mappable<T>& lhs, const GenomicRegion& rhs)
{
    return overlapped_region(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
auto overlapped_region(const ContigRegion& lhs, const Mappable<T>& rhs)
{
    return overlapped_region(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T>
auto overlapped_region(const GenomicRegion& lhs, const Mappable<T>& rhs)
{
    return overlapped_region(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T1, typename T2>
auto overlapped_region(const Mappable<T1>& lhs, const Mappable<T2>& rhs)
{
    return overlapped_region(static_cast<const T1&>(lhs).get_region(),
                             static_cast<const T2&>(rhs).get_region());
}

template <typename T>
auto left_overhang_size(const Mappable<T>& lhs, const ContigRegion& rhs)
{
    return left_overhang_size(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
auto left_overhang_size(const Mappable<T>& lhs, const GenomicRegion& rhs)
{
    return left_overhang_size(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
auto left_overhang_size(const ContigRegion& lhs, const Mappable<T>& rhs)
{
    return left_overhang_size(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T>
auto left_overhang_size(const GenomicRegion& lhs, const Mappable<T>& rhs)
{
    return left_overhang_size(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T1, typename T2>
auto left_overhang_size(const Mappable<T1>& lhs, const Mappable<T2>& rhs)
{
    return left_overhang_size(static_cast<const T1&>(lhs).get_region(),
                              static_cast<const T2&>(rhs).get_region());
}

template <typename T>
auto right_overhang_size(const Mappable<T>& lhs, const ContigRegion& rhs)
{
    return right_overhang_size(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
auto right_overhang_size(const Mappable<T>& lhs, const GenomicRegion& rhs)
{
    return right_overhang_size(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
auto right_overhang_size(const ContigRegion& lhs, const Mappable<T>& rhs)
{
    return right_overhang_size(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T>
auto right_overhang_size(const GenomicRegion& lhs, const Mappable<T>& rhs)
{
    return right_overhang_size(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T1, typename T2>
auto right_overhang_size(const Mappable<T1>& lhs, const Mappable<T2>& rhs)
{
    return right_overhang_size(static_cast<const T1&>(lhs).get_region(),
                               static_cast<const T2&>(rhs).get_region());
}

template <typename T>
auto left_overhang_region(const Mappable<T>& lhs, const ContigRegion& rhs)
{
    return left_overhang_region(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
auto left_overhang_region(const Mappable<T>& lhs, const GenomicRegion& rhs)
{
    return left_overhang_region(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
auto left_overhang_region(const ContigRegion& lhs, const Mappable<T>& rhs)
{
    return left_overhang_region(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T>
auto left_overhang_region(const GenomicRegion& lhs, const Mappable<T>& rhs)
{
    return left_overhang_region(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T1, typename T2>
auto left_overhang_region(const Mappable<T1>& lhs, const Mappable<T2>& rhs)
{
    return left_overhang_region(static_cast<const T1&>(lhs).get_region(),
                                static_cast<const T2&>(rhs).get_region());
}

template <typename T>
auto right_overhang_region(const Mappable<T>& lhs, const ContigRegion& rhs)
{
    return right_overhang_region(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
auto right_overhang_region(const Mappable<T>& lhs, const GenomicRegion& rhs)
{
    return right_overhang_region(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
auto right_overhang_region(const ContigRegion& lhs, const Mappable<T>& rhs)
{
    return right_overhang_region(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T>
auto right_overhang_region(const GenomicRegion& lhs, const Mappable<T>& rhs)
{
    return right_overhang_region(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T1, typename T2>
auto right_overhang_region(const Mappable<T1>& lhs, const Mappable<T2>& rhs)
{
    return right_overhang_region(static_cast<const T1&>(lhs).get_region(),
                                 static_cast<const T2&>(rhs).get_region());
}

template <typename T>
auto closed_region(const Mappable<T>& lhs, const ContigRegion& rhs)
{
    return closed_region(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
auto closed_region(const Mappable<T>& lhs, const GenomicRegion& rhs)
{
    return closed_region(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
auto closed_region(const ContigRegion& lhs, const Mappable<T>& rhs)
{
    return closed_region(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T>
auto closed_region(const GenomicRegion& lhs, const Mappable<T>& rhs)
{
    return closed_region(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T1, typename T2>
auto closed_region(const Mappable<T1>& lhs, const Mappable<T2>& rhs)
{
    return closed_region(static_cast<const T1&>(lhs).get_region(),
                         static_cast<const T2&>(rhs).get_region());
}

template <typename T>
auto head_region(const Mappable<T>& m, const GenomicRegion::SizeType n = 0)
{
    return head_region(static_cast<const T&>(m).get_region(), n);
}

template <typename T>
auto head_position(const Mappable<T>& m)
{
    return head_position(static_cast<const T&>(m).get_region());
}

template <typename T>
auto tail_region(const Mappable<T>& m, const GenomicRegion::SizeType n = 0)
{
    return tail_region(static_cast<const T&>(m).get_region(), n);
}

template <typename T>
auto tail_position(const Mappable<T>& m)
{
    return tail_position(static_cast<const T&>(m).get_region());
}

template <typename T>
auto begin_distance(const Mappable<T>& lhs, const ContigRegion& rhs)
{
    return begin_distance(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
auto begin_distance(const Mappable<T>& lhs, const GenomicRegion& rhs)
{
    return begin_distance(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
auto begin_distance(const ContigRegion& lhs, const Mappable<T>& rhs)
{
    return begin_distance(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T>
auto begin_distance(const GenomicRegion& lhs, const Mappable<T>& rhs)
{
    return begin_distance(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T1, typename T2>
auto begin_distance(const Mappable<T1>& lhs, const Mappable<T2>& rhs)
{
    return begin_distance(static_cast<const T1&>(lhs).get_region(),
                          static_cast<const T2&>(rhs).get_region());
}

template <typename T>
auto end_distance(const Mappable<T>& lhs, const ContigRegion& rhs)
{
    return end_distance(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
auto end_distance(const Mappable<T>& lhs, const GenomicRegion& rhs)
{
    return end_distance(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
auto end_distance(const ContigRegion& lhs, const Mappable<T>& rhs)
{
    return end_distance(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T>
auto end_distance(const GenomicRegion& lhs, const Mappable<T>& rhs)
{
    return end_distance(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T1, typename T2>
auto end_distance(const Mappable<T1>& lhs, const Mappable<T2>& rhs)
{
    return end_distance(static_cast<const T1&>(lhs).get_region(),
                        static_cast<const T2&>(rhs).get_region());
}

#endif
