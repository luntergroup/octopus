//
//  region_utils.h
//  Octopus
//
//  Created by Daniel Cooke on 10/04/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_region_utils_h
#define Octopus_region_utils_h

#include <algorithm> // std::equal_range, std::count_if, std::any_of, std::find_if, std::min_element, std::max_element
#include <cstddef>   // std::size_t
#include <iterator>  // std::distance, std::cbegin, std::cend, std::prev, std::next

#include "genomic_region.h"
#include "mappable.h"

/**
 Returns the leftmost mappable element in the range [first, last).
 
 The range [first, last) is not required to be sorted.
 */
template <typename ForwardIterator>
inline
ForwardIterator leftmost_mappable(ForwardIterator first, ForwardIterator last)
{
    return std::min_element(first, last);
}

/**
 Returns the rightmost mappable element in the range [first, last).
 
 The range [first, last) is not required to be sorted.
 */
template <typename ForwardIterator>
inline
ForwardIterator rightmost_mappable(ForwardIterator first, ForwardIterator last)
{
    return std::max_element(first, last, [] (const auto& lhs, const auto& rhs) {
        return ends_before(lhs, rhs);
    });
}

/**
 Returns the sub-range of Mappable elements in the range [first, last) such that each element
 in the sub-range overlaps a_region.
 
 Requires [first, last) is sorted w.r.t GenomicRegion::operator<
 */
template <typename ForwardIterator>
inline
std::pair<ForwardIterator, ForwardIterator> overlap_range(ForwardIterator first, ForwardIterator last,
                                                          const GenomicRegion& a_region)
{
    return std::equal_range(first, last, a_region,
                            [] (const auto& lhs, const auto& rhs) {
                                return is_before(lhs, rhs);
                            });
}

/**
 Returns the number of mappable elements in the range [first, last) that overlap with a_region.
 
 Requires [first, last) is sorted w.r.t GenomicRegion::operator<
 */
template <typename ForwardIterator>
inline
std::size_t count_overlapped(ForwardIterator first, ForwardIterator last,
                             const GenomicRegion& a_region)
{
    auto overlapped = overlap_range(first, last, a_region);
    return std::distance(overlapped.first, overlapped.second);
}

/**
 Returns the number of Mappable elements in the range [first, last) such that both lhs and rhs overlap
 the the same element.
 
 Requires [first, last) is sorted w.r.t GenomicRegion::operator<
 */
template <typename ForwardIterator>
inline
std::size_t count_shared(ForwardIterator first, ForwardIterator last,
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
inline
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

/**
 Returns the first Mappable element in the range [first2, last2) such that the element overlaps at-least
  one element in the range [first1, last1) with a_region.
 
 Requires [first1, last1) and [first2, last2) are sorted w.r.t GenomicRegion::operator<
 */
template <typename ForwardIterator>
inline
ForwardIterator find_first_shared(ForwardIterator first1, ForwardIterator last1,
                                  ForwardIterator first2, ForwardIterator last2,
                                  const GenomicRegion& a_region)
{
    return std::find_if(first2, last2, [first1, last1, &a_region] (const auto& mappable) {
        return has_shared(first1, last1, mappable, a_region);
    });
}

template <typename ForwardIterator, typename T>
inline
std::pair<ForwardIterator, ForwardIterator> overlap_range(ForwardIterator first, ForwardIterator last,
                                                          const Mappable<T>& m)

{
    return overlap_range(first, last, static_cast<const T&>(m).get_region());
}

template <typename ForwardIterator, typename T>
inline
std::size_t count_overlapped(ForwardIterator first, ForwardIterator last,
                             const Mappable<T>& m)
{
    return count_overlapped(first, last, static_cast<const T&>(m).get_region());
}

template <typename ForwardIterator, typename T1, typename T2>
inline
std::size_t count_shared(ForwardIterator first, ForwardIterator last,
                       const Mappable<T1>& lhs, const Mappable<T2>& rhs)
{
    return count_shared(first, last, static_cast<const T1&>(lhs).get_region(),
                        static_cast<const T2&>(rhs).get_region());
}

template <typename ForwardIterator, typename T1, typename T2>
inline
bool has_shared(ForwardIterator first, ForwardIterator last,
                const Mappable<T1>& lhs, const Mappable<T2>& rhs)
{
    return has_shared(first, last, static_cast<const T1&>(lhs).get_region(),
                      static_cast<const T2&>(rhs).get_region());
}

template <typename ForwardIterator1, typename ForwardIterator2, typename T>
inline
ForwardIterator2 find_first_shared(ForwardIterator1 first1, ForwardIterator1 last1,
                                   ForwardIterator2 first2, ForwardIterator2 last2,
                                   const Mappable<T>& m)
{
    return std::find_if(first2, last2, [first1, last1, &m] (const auto& mappable) {
        return has_shared(first1, last1, mappable, m);
    });
}

/**
 Counts the number of Mappable elements in the range [first2 + 1, last2) that have shared overalapped
 elements in the range [first1, last1) with the first2
 
 Requires [first1, last1) and [first2, last2) to be sorted w.r.t GenomicRegion::operator<
 */
template <typename ForwardIterator1, typename ForwardIterator2>
std::size_t count_if_shared_with_first(ForwardIterator1 first1, ForwardIterator1 last1,
                                       ForwardIterator2 first2, ForwardIterator2 last2)
{
    if (first2 == last2) return 0;
    
    auto first_overlap_range = overlap_range(first1, last1, *first2);
    
    if (first_overlap_range.first == first_overlap_range.second) return 0;
    
    auto last_overlapped_mappable = std::prev(first_overlap_range.second);
    
    auto overlapped_with_last_range = overlap_range(std::next(first2), last2,
                                                    *last_overlapped_mappable);
    
    return std::distance(overlapped_with_last_range.first, overlapped_with_last_range.second);
}

#endif
