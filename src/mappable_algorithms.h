//
//  mappable_algorithms.h
//  Octopus
//
//  Created by Daniel Cooke on 10/04/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_mappable_algorithms_h
#define Octopus_mappable_algorithms_h

#include <algorithm> // std::equal_range, std::is_sorted, std::binary_search, std::count_if, std::any_of,
                     // std::find_if, std::min_element, std::max_element, std::lower_bound, std::find_if_not
                     // std::generate_n, std::transform
#include <numeric>   // std::accumulate
#include <cstddef>   // std::size_t
#include <iterator>  // std::distance, std::cbegin, std::cend, std::prev, std::next, std::reverse_iterator
#include <stdexcept>
#include <boost/iterator/filter_iterator.hpp>
#include <boost/range/iterator_range_core.hpp>

#include "genomic_region.h"
#include "mappable.h"
#include "mappable_ranges.h"

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
 Returns the mappable element in the range [first, last) with the largest size.
 
 The range [first, last) is not required to be sorted.
 */
template <typename ForwardIterator>
inline
ForwardIterator largest_element(ForwardIterator first, ForwardIterator last)
{
    return std::max_element(first, last, [] (const auto& lhs, const auto& rhs) {
        return size(lhs) < size(rhs);
    });
}

/**
 Returns the mappable element in the range [first, last) with the smallest size.
 
 The range [first, last) is not required to be sorted.
 */
template <typename ForwardIterator>
inline
ForwardIterator smallest_element(ForwardIterator first, ForwardIterator last)
{
    return std::min_element(first, last, [] (const auto& lhs, const auto& rhs) {
        return size(lhs) < size(rhs);
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
    return std::lower_bound(first, last, next_position(mappable));
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
std::vector<boost::iterator_range<ForwardIterator>>
bidirectionally_sorted_ranges(ForwardIterator first, ForwardIterator last)
{
    std::vector<boost::iterator_range<ForwardIterator>> result {};
    
    auto it = first;
    
    while (first != last) {
        it = is_bidirectionally_sorted_until(first, last);
        result.emplace_back(first, it);
        first = it;
    }
    
    result.shrink_to_fit();
    
    return result;
}

/**
 Returns the sub-range(s) of Mappable elements in the range [first, last) such that each element
 in the sub-range overlaps mappable.
 
 The returned OverlapRange is a range of filter iterators (i.e. skips over non-overlapped elements).
 
 The algorithm takes linear time if is_bidirectional=false and logorithmic time if is_bidirectional=true.
 
 Requires [first, last) is sorted w.r.t GenomicRegion::operator<
 */
template <typename ForwardIterator, typename MappableType>
inline
OverlapRange<ForwardIterator>
overlap_range(ForwardIterator first, ForwardIterator last, const MappableType& mappable,
              MappableRangeOrder order=MappableRangeOrder::ForwardSorted)
{
    if (order == MappableRangeOrder::BidirectionallySorted) {
        auto range = std::equal_range(first, last, mappable,
                                      [] (const auto& lhs, const auto& rhs) {
                                          return is_before(lhs, rhs);
                                      });
        
        return make_overlap_range(range.first, range.second, mappable);
    }
    
    auto it = find_first_after(first, last, mappable);
    
    return make_overlap_range(std::find_if(first, it, [&mappable] (const auto& m) {
                                            return overlaps(m, mappable); }), it, mappable);
}

/**
 Returns the sub-range(s) of Mappable elements in the range [first, last) such that each element
 in the sub-range overlaps a_region.
 
 The returned OverlapRange is a range of filter iterators (i.e. skips over non-overlapped elements).
 
 Requires [first, last) is sorted w.r.t GenomicRegion::operator<
 */
template <typename ForwardIterator, typename MappableType>
inline
OverlapRange<ForwardIterator>
overlap_range(ForwardIterator first, ForwardIterator last, const MappableType& mappable,
              GenomicRegion::SizeType max_mappable_size)
{
    using MappableType2 = typename ForwardIterator::value_type;
    
    auto it = find_first_after(first, last, mappable);
    
    auto it2 = std::lower_bound(first, it, shift(get_region(mappable), -std::min(get_begin(mappable), max_mappable_size)),
                                [] (const auto& lhs, const auto& rhs) {
                                    return begins_before(lhs, rhs);
                                });
    
    it2 = std::find_if(it2, it, [&mappable] (const auto& m) { return overlaps(m, mappable); });
    
    return make_overlap_range(it2, it, mappable);
}

/**
 Returns true if any of the mappable elements in the range [first, last) overlap with mappable.
 
 Requires [first, last) is sorted w.r.t GenomicRegion::operator<
 */
template <typename BidirectionalIterator, typename MappableType>
inline
bool has_overlapped(BidirectionalIterator first, BidirectionalIterator last,
                    const MappableType& mappable, MappableRangeOrder order=MappableRangeOrder::ForwardSorted)
{
    if (order == MappableRangeOrder::BidirectionallySorted) {
        return std::binary_search(first, last, mappable,
                                  [] (const auto& lhs, const auto& rhs) {
                                      return is_before(lhs, rhs);
                                  });
    } else {
        auto it = find_first_after(first, last, mappable);
        
        using ReverseIterator = std::reverse_iterator<BidirectionalIterator>;
        
        // searches in reverse order on the assumption regions closer to the boundry with
        // mappable are more likely to overlap with mappable.
        
        return std::any_of(ReverseIterator(it), ReverseIterator(std::next(first)),
                           [&mappable] (const auto& m) {
                               return overlaps(mappable, m);
                           }) || overlaps(*first, mappable);
    }
}

/**
 Returns the number of mappable elements in the range [first, last) that overlap with mappable.
 
 Requires [first, last) is sorted w.r.t GenomicRegion::operator<
 */
template <typename ForwardIterator, typename MappableType>
inline
std::size_t count_overlapped(ForwardIterator first, ForwardIterator last, const MappableType& mappable,
                             MappableRangeOrder order=MappableRangeOrder::ForwardSorted)
{
    auto overlapped = overlap_range(first, last, mappable, order);
    return std::distance(overlapped.begin(), overlapped.end());
}

/**
 Returns the number of mappable elements in the range [first, last) that overlap with mappable.
 
 Requires [first, last) is sorted w.r.t GenomicRegion::operator<
 */
template <typename ForwardIterator, typename MappableType>
inline
std::size_t count_overlapped(ForwardIterator first, ForwardIterator last, const MappableType& mappable,
                             GenomicRegion::SizeType max_mappable_size)
{
    auto overlapped = overlap_range(first, last, mappable, max_mappable_size);
    return std::distance(overlapped.begin(), overlapped.end());
}

/**
 Returns the sub-range of Mappable elements in the range [first, last) such that each element
 in the sub-range is contained within mappable.
 
 Requires [first, last) is sorted w.r.t GenomicRegion::operator<
 */
template <typename BidirectionalIterator, typename MappableType>
inline
ContainedRange<BidirectionalIterator>
contained_range(BidirectionalIterator first, BidirectionalIterator last, const MappableType& mappable)
{
    auto it = std::lower_bound(first, last, mappable,
                               [] (const auto& lhs, const auto& rhs) {
                                    return begins_before(lhs, rhs);
                                });
    
    using ReverseIterator = std::reverse_iterator<BidirectionalIterator>;
    
    auto rit = std::find_if(ReverseIterator(find_first_after(it, last, mappable)), ReverseIterator(std::next(it)),
                            [&mappable] (const auto& m) { return contains(mappable, m); });
    
    return make_contained_range(it, rit.base(), mappable);
}

/**
 Returns true if any of the mappable elements in the range [first, last) are contained within mappable.
 
 Requires [first, last) is sorted w.r.t GenomicRegion::operator<
 */
template <typename BidirectionalIterator, typename MappableType>
inline
bool has_contained(BidirectionalIterator first, BidirectionalIterator last, const MappableType& mappable)
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
template <typename BidirectionalIterator, typename MappableType>
inline
std::size_t count_contained(BidirectionalIterator first, BidirectionalIterator last, const MappableType& mappable)
{
    auto contained = contained_range(first, last, mappable);
    return std::distance(contained.begin(), contained.end());
}

/**
 Returns the number of Mappable elements in the range [first, last) that both lhs and rhs overlap.
 
 Requires [first, last) is sorted w.r.t GenomicRegion::operator<
 */
template <typename ForwardIterator, typename MappableType1, typename MappableType2>
inline
std::size_t count_shared(ForwardIterator first, ForwardIterator last,
                         const MappableType1& lhs, const MappableType2& rhs,
                         MappableRangeOrder order=MappableRangeOrder::ForwardSorted)
{
    auto lhs_overlapped = overlap_range(first, last, lhs, order);
    auto rhs_overlapped = overlap_range(first, last, rhs, order);
    
    return (size(lhs_overlapped) <= size(rhs_overlapped)) ?
            std::count_if(lhs_overlapped.begin(), lhs_overlapped.end(),
                          [&rhs] (const auto& region) {
                              return overlaps(region, rhs);
                          }) :
            std::count_if(rhs_overlapped.begin(), rhs_overlapped.end(),
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
                const MappableType1& lhs, const MappableType2& rhs,
                MappableRangeOrder order=MappableRangeOrder::ForwardSorted)
{
    auto lhs_overlapped = overlap_range(first, last, lhs, order);
    auto rhs_overlapped = overlap_range(first, last, rhs, order);
    
    return (size(lhs_overlapped) <= size(rhs_overlapped)) ?
            std::any_of(lhs_overlapped.begin(), lhs_overlapped.end(),
                        [&rhs] (const auto& region) {
                            return overlaps(region, rhs);
                        }) :
            std::any_of(rhs_overlapped.begin(), rhs_overlapped.end(),
                        [&lhs] (const auto& region) {
                            return overlaps(region, lhs);
                        });
}

/**
 Returns the first Mappable element in the range [first2, last2) such that the element shares a region 
 in the range [first1, last1) with mappable.
 
 Requires [first1, last1) and [first2, last2) are sorted w.r.t GenomicRegion::operator<
 */
template <typename ForwardIterator1, typename ForwardIterator2, typename MappableType>
inline
ForwardIterator2 find_first_shared(ForwardIterator1 first1, ForwardIterator1 last1,
                                   ForwardIterator2 first2, ForwardIterator2 last2,
                                   const MappableType& mappable,
                                   MappableRangeOrder order=MappableRangeOrder::ForwardSorted)
{
    return std::find_if(first2, last2,
                        [first1, last1, &mappable, order] (const auto& m) {
                            return has_shared(first1, last1, m, mappable, order);
                        });
}

/**
 Counts the number of Mappable elements in the range [first2 + 1, last2) that have shared overalapped
 elements in the range [first1, last1) with the first2
 
 Requires [first1, last1) and [first2, last2) to be sorted w.r.t GenomicRegion::operator<
 */
template <typename ForwardIterator1, typename ForwardIterator2>
std::size_t count_if_shared_with_first(ForwardIterator1 first1, ForwardIterator1 last1,
                                       ForwardIterator2 first2, ForwardIterator2 last2,
                                       MappableRangeOrder order=MappableRangeOrder::ForwardSorted)
{
    if (first2 == last2) return 0;
    
    auto overlapped = overlap_range(first1, last1, *first2, order);
    
    if (empty(overlapped)) return 0;
    
    return size(overlap_range(std::next(first2), last2, *std::prev(overlapped.end()), order));
}

/**
 Splits a_region into an ordered vector of GenomicRegions of size 1
 */
inline std::vector<GenomicRegion> decompose(const GenomicRegion& a_region)
{
    std::vector<GenomicRegion> result {};
    
    auto num_elements = size(a_region);
    
    if (num_elements == 0) return result;
    
    result.reserve(num_elements);
    
    unsigned n {0};
    
    std::generate_n(std::back_inserter(result), num_elements, [&n, &a_region] () {
        return GenomicRegion {a_region.get_contig_name(), a_region.get_begin() + n, a_region.get_begin() + ++n};
    });
    
    return result;
}

/**
 Returns the GenomicRegion encompassed by the elements in the range [first, last).
 
 Requires [first, last) is sorted w.r.t GenomicRegion::operator<
 */
template <typename ForwardIterator>
inline
GenomicRegion get_encompassing_region(ForwardIterator first, ForwardIterator last)
{
    if (first == last) {
        throw std::runtime_error {"cannot get encompassing region of empty range"};
    }
    
    return get_encompassing(*first, *rightmost_mappable(first, last));
}

/**
 Returns the minimal range of non-overlapping GenomicRegion's such that each element in the range [first, last)
 is contained within a single region.
 
 Requires [first_mappable, last_mappable) is sorted w.r.t GenomicRegion::operator<
 */
template <typename ForwardIterator>
std::vector<GenomicRegion> get_covered_regions(ForwardIterator first, ForwardIterator last)
{
    std::vector<GenomicRegion> result {};
    
    if (first == last) return result;
    
    ForwardIterator first_overlapped {first};
    ForwardIterator rightmost {first};
    
    while (first != last) {
        if (get_begin(*first) > get_end(*rightmost)) {
            result.emplace_back(get_contig_name(*first_overlapped), get_begin(*first_overlapped), get_end(*rightmost));
            rightmost        = first;
            first_overlapped = first;
        } else {
            if (ends_before(*rightmost, *first)) {
                rightmost = first;
            }
        }
        ++first;
    }
    
    result.emplace_back(get_contig_name(*first_overlapped), get_begin(*first_overlapped), get_end(*rightmost));
    
    result.shrink_to_fit();
    
    return result;
}

/**
 Returns all intervening GenomicRegion's between non-overlapping mappables in the range 
 [first_mappable, last_mappable)
 
 Requires [first_mappable, last_mappable) is sorted w.r.t GenomicRegion::operator<
 */
template <typename ForwardIterator>
std::vector<GenomicRegion> get_all_intervening(ForwardIterator first, ForwardIterator last)
{
    if (first == last) return std::vector<GenomicRegion> {};
    
    std::vector<GenomicRegion> result(std::distance(first, last) - 1);
    
    std::transform(first, std::prev(last), std::next(first), result.begin(),
                   [] (const auto& mappable, const auto& next_mappable) {
                       return get_intervening(mappable, next_mappable);
                   });
    
    return result;
}

#endif
