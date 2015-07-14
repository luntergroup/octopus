//
//  region_algorithms.h
//  Octopus
//
//  Created by Daniel Cooke on 10/04/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_region_algorithms_h
#define Octopus_region_algorithms_h

#include <algorithm> // std::equal_range, std::is_sorted, std::binary_search, std::count_if, std::any_of,
                     // std::find_if, std::min_element, std::max_element, std::lower_bound, std::find_if_not
                     // std::generate_n, std::transform
#include <numeric>   // std::accumulate
#include <cstddef>   // std::size_t
#include <iterator>  // std::distance, std::cbegin, std::cend, std::prev, std::next
#include <stdexcept>
#include <boost/iterator/filter_iterator.hpp>
#include <boost/range/iterator_range_core.hpp>

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
        return (ends_equal(lhs, rhs)) ? begins_before(lhs, rhs) : ends_before(lhs, rhs);
    });
}

/**
 Returns the first Mappable element in the range [first, last) that is_after mappable
 
 Requires [first, last) is sorted w.r.t GenomicRegion::operator<
 */
template <typename ForwardIterator, typename MappableType>
inline
ForwardIterator find_first_after(ForwardIterator first, ForwardIterator last, const MappableType& mappable)
{
    return std::upper_bound(first, last, next_position(mappable));
}

/**
 Returns true if the range of Mappable elements in the range [first, last) is sorted
 w.r.t GenomicRegion::operator<, and satisfies the condition 'if lhs < rhs then end(lhs) < end(rhs)
 */
template <typename ForwardIterator>
inline
bool is_bidirectionally_sorted(ForwardIterator first, ForwardIterator last)
{
    return std::is_sorted(first, last, [] (const auto& lhs, const auto& rhs) {
        return lhs < rhs || ends_before(lhs, rhs);
    });
}

/**
 Returns the an iterator to the first element in the range [first, last) that is not sorted
 according w.r.t GenomicRegion::operator< and 'if lhs < rhs then end(lhs) < end(rhs)
 */
template <typename ForwardIterator>
inline
ForwardIterator is_bidirectionally_sorted_until(ForwardIterator first, ForwardIterator last)
{
    return std::is_sorted_until(first, last, [] (const auto& lhs, const auto& rhs) {
        return lhs < rhs || ends_before(lhs, rhs);
    });
}

/**
 Returns the minimum number of sub-ranges in each within the range [first, last) such that each
 sub-range is bidirectionally_sorted.
 */
template <typename ForwardIterator>
inline
std::vector<std::pair<ForwardIterator, ForwardIterator>>
bidirectionally_sorted_ranges(ForwardIterator first, ForwardIterator last)
{
    std::vector<std::pair<ForwardIterator, ForwardIterator>> result {};
    
    auto it = first;
    
    while (first != last) {
        it = is_bidirectionally_sorted_until(first, last);
        result.emplace_back(std::make_pair(first, it));
        first = it;
    }
    
    result.shrink_to_fit();
    
    return result;
}

template <typename MappableType1, typename MappableType2>
class IsOverlapped
{
public:
    IsOverlapped(const MappableType1& mappable) : mappable_ {mappable} {}
    bool operator()(const MappableType2& m) { return overlaps(mappable_, m); }
private:
    const MappableType1& mappable_;
};

template <typename MappableType1, typename MappableType2, typename Iterator>
using OverlapIterator = boost::filter_iterator<IsOverlapped<MappableType1, MappableType2>, Iterator>;

template <typename MappableType1, typename MappableType2, typename Iterator>
using OverlapRange = boost::iterator_range<OverlapIterator<MappableType1, MappableType2, Iterator>>;

template <typename ForwardIterator, typename MappableType>
inline
OverlapRange<MappableType, typename ForwardIterator::value_type, ForwardIterator>
overlap_range2(ForwardIterator first, ForwardIterator last, const MappableType& mappable)
{
    auto it = find_first_after(first, last, mappable);
    auto f = std::find_if(first, it, [&mappable] (const auto& m) { return overlaps(m, mappable); });
    return boost::make_iterator_range(boost::make_filter_iterator<IsOverlapped<MappableType, typename ForwardIterator::value_type>>(IsOverlapped<MappableType, typename ForwardIterator::value_type>(mappable), f, it),
                                      boost::make_filter_iterator<IsOverlapped<MappableType, typename ForwardIterator::value_type>>(IsOverlapped<MappableType, typename ForwardIterator::value_type>(mappable), it, it));
}

/**
 Returns the sub-range of Mappable elements in the range [first, last) such that each element
 in the sub-range overlaps a_region.
 
 Requires [first, last) is sorted bidirectionally sorted
 */
template <typename ForwardIterator, typename MappableType>
inline
std::pair<ForwardIterator, ForwardIterator> overlap_range(ForwardIterator first, ForwardIterator last,
                                                          const MappableType& mappable)
{
    return std::equal_range(first, last, mappable,
                            [] (const auto& lhs, const auto& rhs) {
                                return is_before(lhs, rhs);
                            });
}

/**
 Returns true if any of the mappable elements in the range [first, last) overlap with mappable.
 
 Requires [first, last) is sorted bidirectionally sorted
 */
template <typename BidirectionalIterator, typename MappableType>
inline
bool has_overlapped(BidirectionalIterator first, BidirectionalIterator last, const MappableType& mappable)
{
    return std::binary_search(first, last, mappable,
                              [] (const auto& lhs, const auto& rhs) {
                                  return is_before(lhs, rhs);
                              });
}

/**
 Returns the number of mappable elements in the range [first, last) that overlap with mappable.
 
 Requires [first, last) is sorted bidirectionally sorted
 */
template <typename ForwardIterator, typename MappableType>
inline
std::size_t count_overlapped(ForwardIterator first, ForwardIterator last, const MappableType& mappable)
{
    auto overlapped = overlap_range(first, last, mappable);
    return std::distance(overlapped.first, overlapped.second);
}

/**
 Returns all sub-ranges of Mappable elements in the range [first, last) such that each element
 in each of the sub-ranges overlaps mappable, but none the elements between the sub-ranges overlaps with mappable.
 
 Requires [first, last) is sorted w.r.t GenomicRegion::operator<
 */
template <typename ForwardIterator, typename MappableType>
inline
std::vector<std::pair<ForwardIterator, ForwardIterator>> overlap_ranges(ForwardIterator first, ForwardIterator last,
                                                                        const MappableType& mappable)
{
    std::vector<std::pair<ForwardIterator, ForwardIterator>> result {};
    
    auto it = find_first_after(first, last, mappable);
    
    while (first != it) {
        first = std::find_if(first, it, [&mappable] (const auto& m) { return overlaps(m, mappable); });
        last  = std::find_if_not(first, it, [&mappable] (const auto& m) { return overlaps(m, mappable); });
        if (first != last) result.emplace_back(std::make_pair(first, last));
        first = last;
    }
    
    result.shrink_to_fit();
    
    return result;
}

/**
 Returns true if any of the mappable elements in the range [first, last) overlap with mappable.
 
 Requires [first, last) is sorted w.r.t GenomicRegion::operator<
 */
template <typename BidirectionalIterator, typename MappableType>
inline
bool has_any_overlapped(BidirectionalIterator first, BidirectionalIterator last, const MappableType& mappable)
{
    auto it = find_first_after(first, last, mappable);
    
    using ReverseIterator = std::reverse_iterator<BidirectionalIterator>;
    
    // searches in reverse order on the assumption regions closer to the boundry with
    // mappable are more likely to overlap with mappable. This assumption may not hold if
    // the elements in the range [first, last) have highly variant sizes.
    
    return std::any_of(ReverseIterator(it), ReverseIterator(first++),
                       [&mappable] (const auto& m) {
                           return overlaps(mappable, m);
                       }) || overlaps(*first, mappable);
}

/**
 Returns the number of mappable elements in the range [first, last) that overlap with mappable.
 
 Requires [first, last) is sorted w.r.t GenomicRegion::operator<
 */
template <typename ForwardIterator, typename MappableType>
inline
std::size_t count_all_overlapped(ForwardIterator first, ForwardIterator last, const MappableType& mappable)
{
    auto overlapped = overlap_ranges(first, last, mappable);
    return std::accumulate(std::cbegin(overlapped), std::cend(overlapped), 0,
                           [] (std::size_t lhs, const auto& rhs) {
                               return lhs + std::distance(rhs.first, rhs.second);
                           });
}

/**
 Returns the sub-range of Mappable elements in the range [first, last) such that each element
 in the sub-range is contained within mappable.
 
 Requires [first, last) is sorted w.r.t GenomicRegion::operator<
 */
template <typename ForwardIterator, typename MappableType>
inline
std::pair<ForwardIterator, ForwardIterator> contained_range(ForwardIterator first, ForwardIterator last,
                                                            const MappableType& mappable)
{
    auto it = std::lower_bound(first, last, mappable,
                               [] (const auto& lhs, const auto& rhs) {
                                    return begins_before(lhs, rhs);
                                });
    
    return std::make_pair(it, std::find_if_not(it, last,
                                               [&mappable] (const auto& m) {
                                                   return get_end(m) <= get_end(mappable);
                                               }));
}

/**
 Returns true if any of the mappable elements in the range [first, last) are contained within mappable.
 
 Requires [first, last) is sorted w.r.t GenomicRegion::operator<
 */
template <typename ForwardIterator, typename MappableType>
inline
bool has_contained(ForwardIterator first, ForwardIterator last, const MappableType& mappable)
{
    auto it = std::lower_bound(first, last, mappable,
                               [] (const auto& lhs, const auto& rhs) {
                                    return begins_before(lhs, rhs);
                                });
    
    return (it != last) && get_end(*it) <= get_end(mappable);
}

/**
 Returns the number of mappable elements in the range [first, last) that are contained within mappable.
 
 Requires [first, last) is sorted w.r.t GenomicRegion::operator<
 */
template <typename ForwardIterator, typename MappableType>
inline
std::size_t count_contained(ForwardIterator first, ForwardIterator last, const MappableType& mappable)
{
    auto contained = contained_range(first, last, mappable);
    return std::distance(contained.first, contained.second);
}

/**
 Returns the number of Mappable elements in the range [first, last) such that both lhs and rhs overlap
 the the same element.
 
 Requires [first, last) is sorted w.r.t GenomicRegion::operator<
 */
template <typename ForwardIterator, typename MappableType1, typename MappableType2>
inline
std::size_t count_shared(ForwardIterator first, ForwardIterator last,
                         const MappableType1& lhs, const MappableType2& rhs)
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
template <typename ForwardIterator, typename MappableType1, typename MappableType2>
inline
bool has_shared(ForwardIterator first, ForwardIterator last,
                const MappableType1& lhs, const MappableType2& rhs)
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
template <typename ForwardIterator1, typename ForwardIterator2, typename MappableType>
inline
ForwardIterator2 find_first_shared(ForwardIterator1 first1, ForwardIterator1 last1,
                                  ForwardIterator2 first2, ForwardIterator2 last2,
                                  const MappableType& mappable)
{
    return std::find_if(first2, last2, [first1, last1, &mappable] (const auto& m) {
        return has_shared(first1, last1, m, mappable);
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
    
    auto overlapped_with_last_range = overlap_range(std::next(first2), last2, *std::prev(first_overlap_range.second));
    
    return std::distance(overlapped_with_last_range.first, overlapped_with_last_range.second);
}

/**
 Returns the GenomicRegion encompassed by the elements in the range [first, last).
 
 Requires [first, last) is sorted w.r.t GenomicRegion::operator<
 */
template <typename ForwardIterator>
inline
GenomicRegion encompassing(ForwardIterator first, ForwardIterator last)
{
    if (std::distance(first, last) == 0) {
        throw std::runtime_error {"cannot get encompassed region of empty range"};
    }
    
    return get_encompassing(*first, *rightmost_mappable(first, last));
}

/**
 Splits a_region into an ordered vector of GenomicRegions of size 1
 */
inline std::vector<GenomicRegion> decompose(const GenomicRegion& a_region)
{
    std::vector<GenomicRegion> result {};
    
    auto num_elements = size(a_region);
    
    result.reserve(num_elements);
    
    unsigned n {0};
    
    std::generate_n(std::back_inserter(result), num_elements, [&n, &a_region] () {
        return GenomicRegion {a_region.get_contig_name(), a_region.get_begin() + n, a_region.get_begin() + ++n};
    });
    
    return result;
}

/**
 Returns all overlapped GenomicRegion's in the range [first_mappable, last_mappable) such that
 the returned regions covers the same overall area as the input range (i.e. the maximal overlapping
 regions)
 
 Requires [first_mappable, last_mappable) is sorted w.r.t GenomicRegion::operator<
 */
template <typename ForwardIterator>
std::vector<GenomicRegion> get_all_overlapped(ForwardIterator first_mappable, ForwardIterator last_mappable)
{
    if (first_mappable == last_mappable) return std::vector<GenomicRegion> {};
    
    std::vector<GenomicRegion> result {};
    
    ForwardIterator last_overlapped;
    
    while (first_mappable != last_mappable) {
        last_overlapped = overlap_range(first_mappable, last_mappable, *first_mappable).second;
        result.emplace_back(encompassing(first_mappable, last_overlapped));
        first_mappable = last_overlapped;
    }
    
    return result;
}

/**
 Returns all intervening GenomicRegion's between non-overlapping mappables in the range 
 [first_mappable, last_mappable)
 
 Requires [first_mappable, last_mappable) is sorted w.r.t GenomicRegion::operator<
 */
template <typename ForwardIterator>
std::vector<GenomicRegion> get_all_intervening(ForwardIterator first_mappable, ForwardIterator last_mappable)
{
    if (first_mappable == last_mappable) return std::vector<GenomicRegion> {};
    
    std::vector<GenomicRegion> result(std::distance(first_mappable, last_mappable) - 1);
    
    std::transform(first_mappable, std::prev(last_mappable), std::next(first_mappable), result.begin(),
                   [] (const auto& mappable, const auto& next_mappable) {
                       return get_intervening(mappable, next_mappable);
                   });
    
    return result;
}

#endif
