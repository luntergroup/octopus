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
inline bool begins_before(const Mappable<T>& lhs, const GenomicRegion& rhs)
{
    return begins_before(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
inline bool begins_before(const GenomicRegion& lhs, const Mappable<T>& rhs)
{
    return begins_before(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T>
inline bool begins_before(const Mappable<T>& lhs, const Mappable<T>& rhs)
{
    return begins_before(static_cast<const T&>(lhs).get_region(), static_cast<const T&>(rhs).get_region());
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

template <typename T>
inline bool ends_before(const Mappable<T>& lhs, const Mappable<T>& rhs)
{
    return ends_before(static_cast<const T&>(lhs).get_region(), static_cast<const T&>(rhs).get_region());
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

template <typename T>
inline bool is_before(const Mappable<T>& lhs, const Mappable<T>& rhs) noexcept
{
    return is_before(static_cast<const T&>(lhs).get_region(), static_cast<const T&>(rhs).get_region());
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

template <typename T>
inline bool is_after(const Mappable<T>& lhs, const Mappable<T>& rhs) noexcept
{
    return is_after(static_cast<const T&>(lhs).get_region(), static_cast<const T&>(rhs).get_region());
}

template <typename T>
inline GenomicRegion::DifferenceType overlap_size(const Mappable<T>& lhs, const GenomicRegion& rhs) noexcept
{
    return overlaps(static_cast<const T&>(lhs).get_region(), rhs);
}

template <typename T>
inline GenomicRegion::DifferenceType overlap_size(const GenomicRegion& lhs, const Mappable<T>& rhs) noexcept
{
    return overlaps(lhs, static_cast<const T&>(rhs).get_region());
}

template <typename T>
inline GenomicRegion::DifferenceType overlap_size(const Mappable<T>& lhs, const Mappable<T>& rhs) noexcept
{
    return overlaps(static_cast<const T&>(lhs).get_region(), static_cast<const T&>(rhs).get_region());
}

template <typename T>
inline bool overlaps(const Mappable<T>& lhs, const GenomicRegion& rhs) noexcept
{
    return overlap_size(static_cast<const T&>(lhs), rhs) > 0;
}

template <typename T>
inline bool overlaps(const GenomicRegion& lhs, const Mappable<T>& rhs) noexcept
{
    return overlap_size(lhs, static_cast<const T&>(rhs)) > 0;
}

template <typename T>
inline bool overlaps(const Mappable<T>& lhs, const Mappable<T>& rhs) noexcept
{
    return overlap_size(static_cast<const T&>(lhs), static_cast<const T&>(rhs)) > 0;
}

template <typename T>
inline GenomicRegion::SizeType size(const Mappable<T>& m)
{
    return size(static_cast<const T&>(m).get_region());
}

#endif
