//
//  mappable.h
//  Octopus
//
//  Created by Daniel Cooke on 10/04/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_mappable_h
#define Octopus_mappable_h

#include "genomic_region.h"

/**
 More template black magic. Anything that inherits from Mappable and define the get_region() method 
 will automatically get all the following region comparisons. This is very useful for comparing
 classes that exist on a sequence (i.e. they are mapplable).
 */

template <typename T>
class Mappable {};

template <typename T>
inline bool empty(const Mappable<T>& m)
{
    return empty(static_cast<const T&>(m).get_region());
}

template <typename T>
inline GenomicRegion::SizeType size(const Mappable<T>& m)
{
    return size(static_cast<const T&>(m).get_region());
}

template <typename T>
inline GenomicRegion::StringType get_contig_name(const Mappable<T>& m)
{
    return get_contig_name(static_cast<const T&>(m).get_region());
}

template <typename T>
inline GenomicRegion::SizeType get_begin(const Mappable<T>& m)
{
    return get_begin(static_cast<const T&>(m).get_region());
}

template <typename T>
inline GenomicRegion::SizeType get_end(const Mappable<T>& m)
{
    return get_end(static_cast<const T&>(m).get_region());
}

template <typename T>
inline bool operator==(const Mappable<T>& lhs, const GenomicRegion& rhs)
{
    return static_cast<const T&>(lhs).get_region() == rhs;
}

template <typename T>
inline bool operator==(const GenomicRegion& lhs, const Mappable<T>& rhs)
{
    return lhs == static_cast<const T&>(rhs).get_region();
}

template <typename T>
inline bool operator<(const Mappable<T>& lhs, const GenomicRegion& rhs)
{
    return static_cast<const T&>(lhs).get_region() < rhs;
}

template <typename T>
inline bool operator<(const GenomicRegion& lhs, const Mappable<T>& rhs)
{
    return lhs < static_cast<const T&>(rhs).get_region();
}

template <typename T>
inline bool begins_equal(const GenomicRegion& lhs, const Mappable<T>& rhs)
{
    return begins_equal(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T>
inline bool begins_equal(const Mappable<T>& lhs, const GenomicRegion& rhs)
{
    return begins_equal(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T1, typename T2>
inline bool begins_equal(const Mappable<T1>& lhs, const Mappable<T2>& rhs)
{
    return begins_equal(static_cast<const T1&>(lhs).get_region(), static_cast<const T2&>(rhs).get_region());
}

template <typename T>
inline bool ends_equal(const GenomicRegion& lhs, const Mappable<T>& rhs)
{
    return ends_equal(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T>
inline bool ends_equal(const Mappable<T>& lhs, const GenomicRegion& rhs)
{
    return ends_equal(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T1, typename T2>
inline bool ends_equal(const Mappable<T1>& lhs, const Mappable<T2>& rhs)
{
    return ends_equal(static_cast<const T1&>(lhs).get_region(), static_cast<const T2&>(rhs).get_region());
}

template <typename T>
inline bool begins_before(const Mappable<T>& lhs, const GenomicRegion& rhs)
{
    return begins_before(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
inline bool begins_before(const GenomicRegion& lhs, const Mappable<T>& rhs)
{
    return begins_before(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T1, typename T2>
inline bool begins_before(const Mappable<T1>& lhs, const Mappable<T2>& rhs)
{
    return begins_before(static_cast<const T1&>(lhs).get_region(), static_cast<const T2&>(rhs).get_region());
}

template <typename T>
inline bool ends_before(const Mappable<T>& lhs, const GenomicRegion& rhs)
{
    return ends_before(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
inline bool ends_before(const GenomicRegion& lhs, const Mappable<T>& rhs)
{
    return ends_before(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T1, typename T2>
inline bool ends_before(const Mappable<T1>& lhs, const Mappable<T2>& rhs)
{
    return ends_before(static_cast<const T1&>(lhs).get_region(), static_cast<const T2&>(rhs).get_region());
}

template <typename T>
inline bool is_before(const Mappable<T>& lhs, const GenomicRegion& rhs) noexcept
{
    return is_before(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
inline bool is_before(const GenomicRegion& lhs, const Mappable<T>& rhs) noexcept
{
    return is_before(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T1, typename T2>
inline bool is_before(const Mappable<T1>& lhs, const Mappable<T2>& rhs) noexcept
{
    return is_before(static_cast<const T1&>(lhs).get_region(), static_cast<const T2&>(rhs).get_region());
}

template <typename T>
inline bool is_after(const Mappable<T>& lhs, const GenomicRegion& rhs) noexcept
{
    return is_after(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
inline bool is_after(const GenomicRegion& lhs, const Mappable<T>& rhs) noexcept
{
    return is_after(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T1, typename T2>
inline bool is_after(const Mappable<T1>& lhs, const Mappable<T2>& rhs) noexcept
{
    return is_after(static_cast<const T1&>(lhs).get_region(), static_cast<const T2&>(rhs).get_region());
}

template <typename T>
inline GenomicRegion::DifferenceType overlap_size(const Mappable<T>& lhs, const GenomicRegion& rhs) noexcept
{
    return overlap_size(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
inline GenomicRegion::DifferenceType overlap_size(const GenomicRegion& lhs, const Mappable<T>& rhs) noexcept
{
    return overlap_size(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T1, typename T2>
inline GenomicRegion::DifferenceType overlap_size(const Mappable<T1>& lhs, const Mappable<T2>& rhs) noexcept
{
    return overlap_size(static_cast<const T1&>(lhs).get_region(), static_cast<const T2&>(rhs).get_region());
}

template <typename T>
inline bool overlaps(const Mappable<T>& lhs, const GenomicRegion& rhs) noexcept
{
    return overlaps(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
inline bool overlaps(const GenomicRegion& lhs, const Mappable<T>& rhs) noexcept
{
    return overlaps(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T1, typename T2>
inline bool overlaps(const Mappable<T1>& lhs, const Mappable<T2>& rhs) noexcept
{
    return overlaps(static_cast<const T1&>(lhs).get_region(), static_cast<const T2&>(rhs).get_region());
}

template <typename T>
inline bool contains(const Mappable<T>& lhs, const GenomicRegion& rhs) noexcept
{
    return contains(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
inline bool contains(const GenomicRegion& lhs, const Mappable<T>& rhs) noexcept
{
    return contains(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T1, typename T2>
inline bool contains(const Mappable<T1>& lhs, const Mappable<T2>& rhs) noexcept
{
    return contains(static_cast<const T1&>(lhs).get_region(), static_cast<const T2&>(rhs).get_region());
}

template <typename T>
inline bool are_adjacent(const Mappable<T>& lhs, const GenomicRegion& rhs) noexcept
{
    return are_adjacent(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
inline bool are_adjacent(const GenomicRegion& lhs, const Mappable<T>& rhs) noexcept
{
    return are_adjacent(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T1, typename T2>
inline bool are_adjacent(const Mappable<T1>& lhs, const Mappable<T2>& rhs) noexcept
{
    return are_adjacent(static_cast<const T1&>(lhs).get_region(), static_cast<const T2&>(rhs).get_region());
}

template <typename T>
inline GenomicRegion get_encompassing_region(const Mappable<T>& lhs, const GenomicRegion& rhs)
{
    return get_encompassing_region(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
inline GenomicRegion get_encompassing_region(const GenomicRegion& lhs, const Mappable<T>& rhs)
{
    return get_encompassing_region(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T1, typename T2>
inline GenomicRegion get_encompassing_region(const Mappable<T1>& lhs, const Mappable<T2>& rhs)
{
    return get_encompassing_region(static_cast<const T1&>(lhs).get_region(), static_cast<const T1&>(rhs).get_region());
}

template <typename T>
inline GenomicRegion get_intervening_region(const Mappable<T>& lhs, const GenomicRegion& rhs)
{
    return get_intervening_region(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
inline GenomicRegion get_intervening_region(const GenomicRegion& lhs, const Mappable<T>& rhs)
{
    return get_intervening_region(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T1, typename T2>
inline GenomicRegion get_intervening_region(const Mappable<T1>& lhs, const Mappable<T2>& rhs)
{
    return get_intervening_region(static_cast<const T1&>(lhs).get_region(), static_cast<const T1&>(rhs).get_region());
}

template <typename T>
inline GenomicRegion get_overlapped(const Mappable<T>& lhs, const GenomicRegion& rhs)
{
    return get_overlapped(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
inline GenomicRegion get_overlapped(const GenomicRegion& lhs, const Mappable<T>& rhs)
{
    return get_overlapped(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T1, typename T2>
inline GenomicRegion get_overlapped(const Mappable<T1>& lhs, const Mappable<T2>& rhs)
{
    return get_overlapped(static_cast<const T1&>(lhs).get_region(), static_cast<const T1&>(rhs).get_region());
}

template <typename T>
inline GenomicRegion get_left_overhang(const Mappable<T>& lhs, const GenomicRegion& rhs)
{
    return get_left_overhang(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
inline GenomicRegion get_left_overhang(const GenomicRegion& lhs, const Mappable<T>& rhs)
{
    return get_left_overhang(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T1, typename T2>
inline GenomicRegion get_left_overhang(const Mappable<T1>& lhs, const Mappable<T2>& rhs)
{
    return get_left_overhang(static_cast<const T1&>(lhs).get_region(), static_cast<const T1&>(rhs).get_region());
}

template <typename T>
inline GenomicRegion get_right_overhang(const Mappable<T>& lhs, const GenomicRegion& rhs)
{
    return get_right_overhang(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
inline GenomicRegion get_right_overhang(const GenomicRegion& lhs, const Mappable<T>& rhs)
{
    return get_right_overhang(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T1, typename T2>
inline GenomicRegion get_right_overhang(const Mappable<T1>& lhs, const Mappable<T2>& rhs)
{
    return get_right_overhang(static_cast<const T1&>(lhs).get_region(), static_cast<const T1&>(rhs).get_region());
}

#endif
