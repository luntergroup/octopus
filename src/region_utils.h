//
//  region_utils.h
//  Octopus
//
//  Created by Daniel Cooke on 10/04/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_region_utils_h
#define Octopus_region_utils_h

#include <algorithm> // std::equal_range, std::count_if, std::any_of
#include <cstddef>   // std::size_t

#include "genomic_region.h"
#include "mappable.h"

/**
 Returns the sub-range of Mappable elements in the range [first, last) such that each element
 in the sub-range overlaps a_region.
 
 Requires [first, last) is sorted w.r.t GenomicRegion::operator<
 */
template <typename ForwardIterator>
inline std::pair<ForwardIterator, ForwardIterator> overlap_range(ForwardIterator first, ForwardIterator last,
                                                                 const GenomicRegion& a_region)
{
    return std::equal_range(first, last, a_region,
                            [] (const auto& other, const auto& a_region) {
                                return is_before(other, a_region);
                            });
}

/**
 Returns the number of Mappable elements in the range [first, last) such that both lhs and rhs overlap
 the the same region.
 
 Requires [first, last) is sorted w.r.t GenomicRegion::operator<
 */
template <typename ForwardIterator>
std::size_t num_shared(ForwardIterator first, ForwardIterator last,
                       const GenomicRegion& lhs, const GenomicRegion& rhs)
{
    auto lhs_overlap_range = overlap_range(first, last, lhs);
    auto rhs_overlap_range = overlap_range(first, last, lhs);
    
    return (std::distance(lhs_overlap_range.first, lhs_overlap_range.second) <=
            std::distance(rhs_overlap_range.first, rhs_overlap_range.second)) ?
            std::count_if(lhs_overlap_range.first, lhs_overlap_range.second,
                          [&rhs] (const auto& region) {
                              return overlaps(region, rhs);
                          }) :
            std::count_if(rhs_overlap_range.first, rhs_overlap_range.second,
                          [&lhs] (const auto& region) {
                              return overlaps(region, lhs);
                          });
}

/**
 Returns if any of the Mappable elements in the range [first, last) overlaps both lhs and rhs.
 
 Requires [first, last) is sorted w.r.t GenomicRegion::operator<
 */
template <typename ForwardIterator>
bool has_shared(ForwardIterator first, ForwardIterator last,
                const GenomicRegion& lhs, const GenomicRegion& rhs)
{
    auto lhs_overlap_range = overlap_range(first, last, lhs);
    auto rhs_overlap_range = overlap_range(first, last, lhs);
    
    return (std::distance(lhs_overlap_range.first, lhs_overlap_range.second) <=
            std::distance(rhs_overlap_range.first, rhs_overlap_range.second)) ?
            std::any_of(lhs_overlap_range.first, lhs_overlap_range.second,
                        [&rhs] (const auto& region) {
                            return overlaps(region, rhs);
                        }) :
            std::any_of(rhs_overlap_range.first, rhs_overlap_range.second,
                        [&lhs] (const auto& region) {
                            return overlaps(region, lhs);
                        });
}

template <typename ForwardIterator, typename T>
inline std::pair<ForwardIterator, ForwardIterator> overlap_range(ForwardIterator first, ForwardIterator last,
                                                                 const Mappable<T>& m)

{
    return overlap_range(first, last, static_cast<const T&>(m).get_region());
}

template <typename ForwardIterator, typename T>
std::size_t num_shared(ForwardIterator first, ForwardIterator last,
                       const Mappable<T>& lhs, const Mappable<T>& rhs)
{
    return num_shared(first, last, static_cast<const T&>(lhs).get_region(), static_cast<const T&>(rhs).get_region());
}

template <typename ForwardIterator, typename T>
bool has_shared(ForwardIterator first, ForwardIterator last,
                const GenomicRegion& lhs, const GenomicRegion& rhs)
{
    return has_shared(first, last, static_cast<const T&>(lhs).get_region(), static_cast<const T&>(rhs).get_region());
}

#endif
