//
//  region_utils.h
//  Octopus
//
//  Created by Daniel Cooke on 10/04/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_region_utils_h
#define Octopus_region_utils_h

#include <algorithm> // std::equal_range, std::count_if, std::any_of, std::find_if, std::transform, std::min_element
#include <cstddef>   // std::size_t
#include <iterator>  // std::distance, std::cbegin etc
#include <vector>

#include "genomic_region.h"
#include "mappable.h"

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
                            [] (const auto& other, const auto& a_region) {
                                return is_before(other, a_region);
                            });
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

template <typename MappableMap, typename Mappable>
inline
GenomicRegion leftmost_overlapping(const MappableMap& mappables, const Mappable& mappable)
{
    auto result = mappable.get_region();
    using Iterator = typename MappableMap::mapped_type::const_iterator;
    std::pair<Iterator, Iterator> current;
    
    for (const auto& map_pair : mappables) {
        current = overlap_range(std::cbegin(map_pair.second), std::cend(map_pair.second), mappable);
        
        if (current.first != std::cend(map_pair.second) && current.first->get_region() < result) {
            result = current.first->get_region();
        }
    }
    
    return result;
}

template <typename MappableMap, typename Mappable>
inline
GenomicRegion rightmost_overlapping(const MappableMap& mappables, const Mappable& mappable)
{
    auto result = mappable.get_region();
    using Iterator = typename MappableMap::mapped_type::const_iterator;
    std::pair<Iterator, Iterator> current;
    
    for (const auto& map_pair : mappables) {
        current = overlap_range(std::cbegin(map_pair.second), std::cend(map_pair.second), mappable);
        
        if (current.first != current.second && std::prev(current.second)->get_region() > result) {
            result = std::prev(current.second)->get_region();
        }
    }
    
    return result;
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

template <typename MappableMap, typename Mappable1, typename Mappable2>
inline
bool has_shared(const MappableMap& mappables, const Mappable1& lhs, const Mappable2& rhs)
{
    return std::any_of(std::cbegin(mappables), std::cend(mappables),
                       [&lhs, &rhs] (const auto& map_pair) {
                           return has_shared(std::cbegin(map_pair.second), std::cend(map_pair.second),
                                             lhs, rhs);
                       });
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
 Works the same way as the other find_first_shared, but operators over a map of ranges
 (e.g. std::map or std::unordered_map).
 */
template <typename ForwardIterator, typename MappableMap, typename Mappable>
ForwardIterator find_first_shared(const MappableMap& mappables,
                                  ForwardIterator first, ForwardIterator last,
                                  const Mappable& mappable)
{
    if (mappables.empty()) return last;
    
    if (mappables.size() == 1) {
        return find_first_shared(std::cbegin(std::cbegin(mappables)->second), std::cend(std::cbegin(mappables)->second),
                                 first, last, mappable);
    }
    
    std::vector<ForwardIterator> smallest {};
    smallest.reserve(mappables.size());
    
    std::transform(std::cbegin(mappables), std::cend(mappables), std::back_inserter(smallest),
                   [first, last, &mappable] (const auto& p) {
                       return find_first_shared(std::cbegin(p.second), std::cend(p.second), first, last,
                                                mappable);
                   });
    
    return *std::min_element(std::cbegin(smallest), std::cend(smallest), []
                             (ForwardIterator lhs, ForwardIterator rhs) {
                                 return *lhs < *rhs;
                             });
}

/**
 Counts the number of Mappable elements in the range [first2 + 1, last2) that have shared overalapped
 elements in the range [first1, last1) with the first2
 */
template <typename ForwardIterator1, typename ForwardIterator2>
std::size_t count_if_shared_with_first(ForwardIterator1 first1, ForwardIterator1 last1,
                                       ForwardIterator2 first2, ForwardIterator2 last2)
{
    if (first2 == last2) return 0;
    
    auto first_overlap_range = overlap_range(first1, last1, *first2);
    
    if (first_overlap_range.first == first_overlap_range.second) return 0;
    
    auto last_overlapped_mappable = std::prev(first_overlap_range.second);
    
    auto overlapped_with_last_range = overlap_range(std::next(first2), last2, *last_overlapped_mappable);
    
    return std::distance(overlapped_with_last_range.first, overlapped_with_last_range.second);
}

template <typename MappableMap, typename ForwardIterator>
inline
std::size_t max_count_if_shared_with_first(const MappableMap& mappables,
                                           ForwardIterator first, ForwardIterator last)
{
    std::size_t maximum {0};
    std::size_t count {};
    
    for (const auto& map_pair : mappables) {
        count = count_if_shared_with_first(std::cbegin(map_pair.second), std::cend(map_pair.second),
                                           first, last);
        if (count > maximum) {
            maximum = count;
        }
    }
    
    return maximum;
}

#endif
