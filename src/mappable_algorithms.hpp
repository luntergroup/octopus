//
//  mappable_algorithms.hpp
//  Octopus
//
//  Created by Daniel Cooke on 10/04/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_mappable_algorithms_hpp
#define Octopus_mappable_algorithms_hpp

#include <algorithm> // std::equal_range, std::is_sorted, std::binary_search, std::count_if, std::any_of,
                     // std::find_if, std::min_element, std::max_element, std::lower_bound, std::find_if_not
                     // std::generate_n, std::transform
#include <numeric>   // std::accumulate
#include <cstddef>   // size_t
#include <iterator>  // std::distance, std::cbegin, std::cend, std::prev, std::next,
                     // std::make_reverse_iterator, std::iterator_traits
#include <stdexcept>
#include <boost/iterator/filter_iterator.hpp>
#include <boost/range/iterator_range_core.hpp>

#include "contig_region.hpp"
#include "genomic_region.hpp"
#include "mappable.hpp"
#include "mappable_ranges.hpp"

#include <iostream> // DEBUG

template <typename ForwardIterator>
GenomicRegion::SizeType sum_sizes(ForwardIterator first, ForwardIterator last)
{
    return std::accumulate(first, last, GenomicRegion::SizeType {},
                           [] (const auto curr, const auto& mappable) { return curr + size(mappable); });
}

template <typename Container>
GenomicRegion::SizeType sum_sizes(const Container& mappables)
{
    return sum_sizes(std::cbegin(mappables), std::cend(mappables));
}

/**
 Returns the leftmost mappable element in the range [first, last).
 
 The range [first, last) is not required to be sorted.
 */
template <typename ForwardIterator>
ForwardIterator leftmost_mappable(ForwardIterator first, ForwardIterator last)
{
    return std::min_element(first, last);
}

template <typename Container>
typename Container::const_iterator leftmost_mappable(const Container& mappables)
{
    return leftmost_mappable(std::cbegin(mappables), std::cend(mappables));
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

template <typename Container>
typename Container::const_iterator rightmost_mappable(const Container& mappables)
{
    return rightmost_mappable(std::cbegin(mappables), std::cend(mappables));
}

/**
 Returns the mappable element in the range [first, last) with the largest size.
 
 The range [first, last) is not required to be sorted.
 */
template <typename ForwardIterator>
inline
ForwardIterator largest_element(ForwardIterator first, ForwardIterator last)
{
    return std::max_element(first, last,
                            [] (const auto& lhs, const auto& rhs) { return size(lhs) < size(rhs); });
}

template <typename Container>
typename Container::const_iterator largest_element(const Container& mappables)
{
    return largest_element(std::cbegin(mappables), std::cend(mappables));
}

/**
 Returns the mappable element in the range [first, last) with the smallest size.
 
 The range [first, last) is not required to be sorted.
 */
template <typename ForwardIterator>
inline
ForwardIterator smallest_element(ForwardIterator first, ForwardIterator last)
{
    return std::min_element(first, last,
                            [] (const auto& lhs, const auto& rhs) { return size(lhs) < size(rhs); });
}

template <typename Container>
typename Container::const_iterator smallest_element(const Container& mappables)
{
    return smallest_element(std::cbegin(mappables), std::cend(mappables));
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

template <typename Container, typename MappableType>
typename Container::const_iterator find_first_after(const Container& mappables, const MappableType& mappable)
{
    return find_first_after(std::cbegin(mappables), std::cend(mappables), mappable);
}

/**
 Returns true if the range of Mappable elements in the range [first, last) is sorted
 w.r.t GenomicRegion::operator<, and satisfies the condition 'if lhs < rhs then end(lhs) < end(rhs)
 */
template <typename ForwardIterator>
inline
bool is_bidirectionally_sorted(ForwardIterator first, ForwardIterator last)
{
    return std::is_sorted(first, last,
                          [] (const auto& lhs, const auto& rhs) { return lhs < rhs || ends_before(lhs, rhs); });
}

template <typename Container>
bool is_bidirectionally_sorted(const Container& mappables)
{
    return is_bidirectionally_sorted(std::cbegin(mappables), std::cend(mappables));
}

/**
 Returns the an iterator to the first element in the range [first, last) that is not sorted
 according w.r.t GenomicRegion::operator< and 'if lhs < rhs then end(lhs) < end(rhs)
 */
template <typename ForwardIterator>
inline
ForwardIterator is_bidirectionally_sorted_until(ForwardIterator first, ForwardIterator last)
{
    return std::is_sorted_until(first, last,
                                [] (const auto& lhs, const auto& rhs) { return lhs < rhs || ends_before(lhs, rhs); });
}

template <typename Container>
typename Container::const_iterator is_bidirectionally_sorted_until(const Container& mappables)
{
    return is_bidirectionally_sorted_until(std::cbegin(mappables), std::cend(mappables));
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
OverlapRange<ForwardIterator>
overlap_range(ForwardIterator first, ForwardIterator last, const MappableType& mappable,
              MappableRangeOrder order = MappableRangeOrder::ForwardSorted)
{
    if (order == MappableRangeOrder::BidirectionallySorted) {
        auto overlapped = std::equal_range(first, last, mappable,
                                           [] (const auto& lhs, const auto& rhs) {
                                               return is_before(lhs, rhs);
                                           });
        
        // now we need to try and push these boundries out as the range does not fully capture
        // insertions
        
        overlapped.second = std::find_if_not(overlapped.second, last,
                                             [&mappable] (const auto& m) { return overlaps(m, mappable); });
        
        auto it = std::find_if_not(std::make_reverse_iterator(overlapped.first),
                                   std::make_reverse_iterator(first),
                                   [&mappable] (const auto& m) { return overlaps(m, mappable); });
        
        return make_overlap_range(it.base(), overlapped.second, mappable);
    }
    
    const auto it = find_first_after(first, last, mappable);
    
    return make_overlap_range(std::find_if(first, it, [&mappable] (const auto& m) {
                                            return overlaps(m, mappable); }), it, mappable);
}

template <typename Container, typename MappableType>
OverlapRange<typename Container::const_iterator>
overlap_range(const Container& mappables, const MappableType& mappable,
              MappableRangeOrder order = MappableRangeOrder::ForwardSorted)
{
    return overlap_range(std::cbegin(mappables), std::cend(mappables), mappable, order);
}

/**
 Returns the sub-range(s) of Mappable elements in the range [first, last) such that each element
 in the sub-range overlaps a_region.
 
 The returned OverlapRange is a range of filter iterators (i.e. skips over non-overlapped elements).
 
 Requires [first, last) is sorted w.r.t GenomicRegion::operator<
 */
template <typename ForwardIterator, typename MappableType>
OverlapRange<ForwardIterator>
overlap_range(ForwardIterator first, ForwardIterator last, const MappableType& mappable,
              GenomicRegion::SizeType max_mappable_size)
{
    using MappableType2 = typename std::iterator_traits<ForwardIterator>::value_type;
    
    auto it = find_first_after(first, last, mappable);
    
    auto it2 = std::lower_bound(first, it, shift(get_region(mappable),
                                                 -std::min(get_begin(mappable), max_mappable_size)),
                                [] (const auto& lhs, const auto& rhs) { return begins_before(lhs, rhs); });
    
    it2 = std::find_if(it2, it, [&mappable] (const auto& m) { return overlaps(m, mappable); });
    
    return make_overlap_range(it2, it, mappable);
}

template <typename Container, typename MappableType>
OverlapRange<typename Container::const_iterator>
overlap_range(const Container& mappables, const MappableType& mappable,
              GenomicRegion::SizeType max_mappable_size)
{
    return overlap_range(std::cbegin(mappables), std::cend(mappables), mappable, max_mappable_size);
}

template <typename Container, typename MappableType>
Container
copy_overlapped(const Container& mappables, const MappableType& mappable,
                MappableRangeOrder order = MappableRangeOrder::ForwardSorted)
{
    const auto overlapped = overlap_range(mappables, mappable, order);
    return Container {std::cbegin(overlapped), std::cend(overlapped)};
}

template <typename Container, typename MappableType>
Container
copy_overlapped(const Container& mappables, const MappableType& mappable,
                GenomicRegion::SizeType max_mappable_size)
{
    const auto overlapped = overlap_range(mappables, mappable, max_mappable_size);
    return Container {std::cbegin(overlapped), std::cend(overlapped)};
}

/**
 Returns true if any of the mappable elements in the range [first, last) overlap with mappable.
 
 Requires [first, last) is sorted w.r.t GenomicRegion::operator<
 */
template <typename BidirectionalIterator, typename MappableType>
bool has_overlapped(BidirectionalIterator first, BidirectionalIterator last,
                    const MappableType& mappable,
                    MappableRangeOrder order = MappableRangeOrder::ForwardSorted)
{
    if (order == MappableRangeOrder::BidirectionallySorted) {
        return std::binary_search(first, last, mappable,
                                  [] (const auto& lhs, const auto& rhs) { return is_before(lhs, rhs); });
    } else {
        if (first == last) return false;
        
        auto it = find_first_after(first, last, mappable);
        
        // searches in reverse order on the assumption regions closer to the boundry with
        // mappable are more likely to overlap with mappable.
        return std::any_of(std::make_reverse_iterator(it), std::make_reverse_iterator(first),
                                                         [&mappable] (const auto& m) {
                                                             return overlaps(mappable, m);
                                                         });
    }
}

template <typename Container, typename MappableType>
bool has_overlapped(const Container& container, const MappableType& mappable,
                    MappableRangeOrder order = MappableRangeOrder::ForwardSorted)
{
    return has_overlapped(std::cbegin(container), std::cend(container), mappable, order);
}

/**
 Returns the number of mappable elements in the range [first, last) that overlap with mappable.
 
 Requires [first, last) is sorted w.r.t GenomicRegion::operator<
 */
template <typename ForwardIterator, typename MappableType>
size_t count_overlapped(ForwardIterator first, ForwardIterator last, const MappableType& mappable,
                        MappableRangeOrder order = MappableRangeOrder::ForwardSorted)
{
    const auto overlapped = overlap_range(first, last, mappable, order);
    return std::distance(overlapped.begin(), overlapped.end());
}

template <typename Container, typename MappableType>
size_t count_overlapped(const Container& container, const MappableType& mappable,
                        MappableRangeOrder order = MappableRangeOrder::ForwardSorted)
{
    return count_overlapped(std::cbegin(container), std::cend(container), mappable, order);
}

/**
 Returns the number of mappable elements in the range [first, last) that overlap with mappable.
 
 Requires [first, last) is sorted w.r.t GenomicRegion::operator<
 */
template <typename ForwardIterator, typename MappableType>
size_t count_overlapped(ForwardIterator first, ForwardIterator last, const MappableType& mappable,
                        GenomicRegion::SizeType max_mappable_size)
{
    const auto overlapped = overlap_range(first, last, mappable, max_mappable_size);
    return std::distance(overlapped.begin(), overlapped.end());
}

template <typename Container, typename MappableType>
size_t count_overlapped(const Container& container, const MappableType& mappable,
                        GenomicRegion::SizeType max_mappable_size)
{
    return count_overlapped(std::cbegin(container), std::cend(container), mappable, max_mappable_size);
}

/**
 Returns true if the range [first, last) contains the exact region given my mappable
 
 Requires [first, last) is sorted w.r.t GenomicRegion::operator<
 */
template <typename ForwardIterator, typename MappableType>
bool has_exact_overlap(ForwardIterator first, ForwardIterator last, const MappableType& mappable,
                       MappableRangeOrder order = MappableRangeOrder::ForwardSorted)
{
    if (order == MappableRangeOrder::BidirectionallySorted) {
        return std::binary_search(first, last, mappable);
    }
    
    const auto overlapped = overlap_range(first, last, mappable);
    
    return std::find_if(std::cbegin(overlapped), std::cend(overlapped),
                        [&mappable] (const auto& e) { return is_same_region(mappable, e); }) != std::cend(overlapped);
}

/**
 Returns the sub-range of Mappable elements in the range [first, last) such that each element
 in the sub-range is contained within mappable.
 
 Requires [first, last) is sorted w.r.t GenomicRegion::operator<
 */
template <typename BidirectionalIterator, typename MappableType>
ContainedRange<BidirectionalIterator>
contained_range(BidirectionalIterator first, BidirectionalIterator last, const MappableType& mappable)
{
    auto it = std::lower_bound(first, last, mappable,
                               [] (const auto& lhs, const auto& rhs) { return begins_before(lhs, rhs); });
    
    auto it2 = find_first_after(it, last, mappable);
    
    if (it == it2) return make_contained_range(it, it2, mappable);
    
    auto rit = std::find_if(std::make_reverse_iterator(it2), std::make_reverse_iterator(std::next(it)),
                            [&mappable] (const auto& m) { return contains(mappable, m); });
    
    return make_contained_range(it, rit.base(), mappable);
}

template <typename Container, typename MappableType>
ContainedRange<typename Container::const_iterator>
contained_range(const Container& mappables, const MappableType& mappable)
{
    return contained_range(std::cbegin(mappables), std::cend(mappables), mappable);
}

/**
 Returns true if any of the mappable elements in the range [first, last) are contained within mappable.
 
 Requires [first, last) is sorted w.r.t GenomicRegion::operator<
 */
template <typename BidirectionalIterator, typename MappableType>
bool has_contained(BidirectionalIterator first, BidirectionalIterator last, const MappableType& mappable)
{
    auto it = std::lower_bound(first, last, mappable,
                               [] (const auto& lhs, const auto& rhs) { return begins_before(lhs, rhs); });
    
    return (it != last) && get_end(*it) <= get_end(mappable);
}

template <typename Container, typename MappableType>
bool has_contained(const Container& container, const MappableType& mappable)
{
    return has_contained(std::cbegin(container), std::cend(container), mappable);
}

/**
 Returns the number of mappable elements in the range [first, last) that are contained within mappable.
 
 Requires [first, last) is sorted w.r.t GenomicRegion::operator<
 */
template <typename BidirectionalIterator, typename MappableType>
size_t count_contained(BidirectionalIterator first, BidirectionalIterator last, const MappableType& mappable)
{
    const auto contained = contained_range(first, last, mappable);
    return std::distance(contained.begin(), contained.end());
}

template <typename Container, typename MappableType>
size_t count_contained(const Container& container, const MappableType& mappable)
{
    return count_contained(std::cbegin(container), std::cend(container), mappable);
}

/**
 Returns the number of Mappable elements in the range [first, last) that both lhs and rhs overlap.
 
 Requires [first, last) is sorted w.r.t GenomicRegion::operator<
 */
template <typename ForwardIterator, typename MappableType1, typename MappableType2>
size_t count_shared(ForwardIterator first, ForwardIterator last,
                    const MappableType1& lhs, const MappableType2& rhs,
                    MappableRangeOrder order = MappableRangeOrder::ForwardSorted)
{
    const auto lhs_overlapped = overlap_range(first, last, lhs, order);
    const auto rhs_overlapped = overlap_range(first, last, rhs, order);
    
    return (size(lhs_overlapped) <= size(rhs_overlapped)) ?
            std::count_if(lhs_overlapped.begin(), lhs_overlapped.end(),
                          [&rhs] (const auto& region) { return overlaps(region, rhs); }) :
            std::count_if(rhs_overlapped.begin(), rhs_overlapped.end(),
                          [&lhs] (const auto& region) { return overlaps(region, lhs); });
}

template <typename Container, typename MappableType1, typename MappableType2>
size_t count_shared(const Container& container,
                       const MappableType1& lhs, const MappableType2& rhs,
                       MappableRangeOrder order = MappableRangeOrder::ForwardSorted)
{
    return count_shared(std::cbegin(container), std::cend(container), lhs, rhs, order);
}

template <typename ForwardIterator>
bool has_exact(ForwardIterator first, ForwardIterator last, const typename ForwardIterator::value_type& mappable)
{
    const auto contained = contained_range(first, last, mappable);
    return std::find(std::cbegin(contained), std::cend(contained), mappable) == std::cend(contained);
}

template <typename Container>
bool has_exact(const Container& mappables, const typename Container::value_type& mappable)
{
    return has_exact(std::cbegin(mappables), std::cend(mappables), mappable);
}

/**
 Returns if any of the Mappable elements in the range [first, last) overlaps both lhs and rhs.
 
 Requires [first, last) is sorted w.r.t GenomicRegion::operator<
 */
template <typename ForwardIterator, typename MappableType1, typename MappableType2>
bool has_shared(ForwardIterator first, ForwardIterator last,
                const MappableType1& lhs, const MappableType2& rhs,
                MappableRangeOrder order = MappableRangeOrder::ForwardSorted)
{
    const auto lhs_overlapped = overlap_range(first, last, lhs, order);
    const auto rhs_overlapped = overlap_range(first, last, rhs, order);
    
    return (size(lhs_overlapped) <= size(rhs_overlapped)) ?
            std::any_of(lhs_overlapped.begin(), lhs_overlapped.end(),
                        [&rhs] (const auto& region) { return overlaps(region, rhs); }) :
            std::any_of(rhs_overlapped.begin(), rhs_overlapped.end(),
                        [&lhs] (const auto& region) { return overlaps(region, lhs); });
}

/**
 Returns the first Mappable element in the range [first2, last2) such that the element shares a region 
 in the range [first1, last1) with mappable.
 
 Requires [first1, last1) and [first2, last2) are sorted w.r.t GenomicRegion::operator<
 */
template <typename ForwardIterator1, typename ForwardIterator2, typename MappableType>
ForwardIterator2 find_first_shared(ForwardIterator1 first1, ForwardIterator1 last1,
                                   ForwardIterator2 first2, ForwardIterator2 last2,
                                   const MappableType& mappable,
                                   MappableRangeOrder order = MappableRangeOrder::ForwardSorted)
{
    return std::find_if(first2, last2,
                        [first1, last1, &mappable, order] (const auto& m) {
                            return has_shared(first1, last1, m, mappable, order);
                        });
}

/**
 Counts the number of Mappable elements in the range [first2 + 1, last2) that have shared
 elements in the range [first1, last1) with first2
 
 Requires [first1, last1) and [first2, last2) to be sorted w.r.t GenomicRegion::operator<
 */
template <typename ForwardIterator1, typename ForwardIterator2>
size_t count_if_shared_with_first(ForwardIterator1 first1, ForwardIterator1 last1,
                                  ForwardIterator2 first2, ForwardIterator2 last2,
                                  MappableRangeOrder order = MappableRangeOrder::ForwardSorted)
{
    if (first2 == last2) return 0;
    
    const auto overlapped = overlap_range(first1, last1, *first2, order);
    
    if (empty(overlapped)) return 0;
    
    return size(overlap_range(std::next(first2), last2, *std::prev(overlapped.end()), order));
}

template <typename Container>
auto get_regions(const Container& mappables)
{
    std::vector<decltype(get_region(typename Container::value_type()))> result {};
    result.reserve(mappables.size());
    
    for (const auto& mappable : mappables) {
        result.emplace_back(get_region(mappable));
    }
    
    return result;
}

/**
 Splits mappable into an ordered vector of GenomicRegions of size 1
 */
template <typename MappableType>
std::vector<GenomicRegion> decompose(const MappableType& mappable)
{
    std::vector<GenomicRegion> result {};
    
    const auto num_elements = size(mappable);
    
    if (num_elements == 0) return result;
    
    result.reserve(num_elements);
    
    GenomicRegion::SizeType n {0};
    
    std::generate_n(std::back_inserter(result), num_elements, [&n, &mappable] () {
        return GenomicRegion {get_contig_name(mappable), get_begin(mappable) + n, get_begin(mappable) + ++n};
    });
    
    return result;
}

template <typename MappableType>
std::vector<GenomicRegion> decompose(const MappableType& mappable, GenomicRegion::SizeType n)
{
    std::vector<GenomicRegion> result {};
    
    if (n == 0) return result;
    
    const auto num_elements = size(mappable) / n;
    
    if (num_elements == 0) return result;
    
    result.reserve(num_elements);
    
    const auto& contig = get_contig_name(mappable);
    auto curr = get_begin(mappable);
    
    std::generate_n(std::back_inserter(result), num_elements, [&contig, &curr, n] () {
        auto tmp = curr;
        curr += n;
        return GenomicRegion {contig, tmp, tmp + n};
    });
    
    return result;
}

/**
 Returns the GenomicRegion encompassed by the elements in the range [first, last).
 
 Requires [first, last) is sorted w.r.t GenomicRegion::operator<
 */
template <typename ForwardIterator>
GenomicRegion get_encompassing_region(ForwardIterator first, ForwardIterator last)
{
    if (first == last) {
        throw std::runtime_error {"get_encompassing_region given empty range"};
    }
    
    return get_encompassing(*first, *rightmost_mappable(first, last));
}

template <typename Container>
GenomicRegion get_encompassing_region(const Container& mappables)
{
    return get_encompassing_region(std::cbegin(mappables), std::cend(mappables));
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
    
    auto first_overlapped = first;
    auto rightmost        = first;
    
    while (first != last) {
        if (get_begin(*first) > get_end(*rightmost)) {
            result.emplace_back(get_contig_name(*first_overlapped), get_begin(*first_overlapped), get_end(*rightmost));
            rightmost        = first;
            first_overlapped = first;
        } else if (ends_before(*rightmost, *first)) {
            rightmost = first;
        }
        ++first;
    }
    
    result.emplace_back(get_contig_name(*first_overlapped), get_begin(*first_overlapped), get_end(*rightmost));
    
    result.shrink_to_fit();
    
    return result;
}

template <typename Container>
std::vector<GenomicRegion> get_covered_regions(const Container& mappables)
{
    return get_covered_regions(std::cbegin(mappables), std::cend(mappables));
}

/**
 Returns all intervening GenomicRegion's between non-overlapping mappables in the range 
 [first, last)
 
 Requires [first, last) is sorted w.r.t GenomicRegion::operator<
 */
template <typename ForwardIterator>
std::vector<GenomicRegion> get_all_intervening(ForwardIterator first, ForwardIterator last)
{
    if (first == last) return std::vector<GenomicRegion> {};
    
    std::vector<GenomicRegion> result {};
    result.reserve(std::distance(first, last) - 1);
    
    std::transform(first, std::prev(last), std::next(first), std::back_inserter(result),
                   [] (const auto& mappable, const auto& next_mappable) {
                       return get_intervening(mappable, next_mappable);
                   });
    
    return result;
}

template <typename Container>
std::vector<GenomicRegion> get_all_intervening(const Container& mappables)
{
    return get_all_intervening(std::cbegin(mappables), std::cend(mappables));
}

/**
 Returns all intervening GenomicRegion's between non-overlapping mappables in the range
 [first, last), and also any flanking regions of mappable if the range [first, last) is 
 contained within mappable.
 
 Requires [first, last) is sorted w.r.t GenomicRegion::operator<
 */
template <typename ForwardIterator, typename Mappable>
std::vector<GenomicRegion> get_all_intervening(ForwardIterator first, ForwardIterator last, const Mappable& mappable)
{
    if (first == last) return std::vector<GenomicRegion> {};
    
    std::vector<GenomicRegion> result {};
    result.reserve(std::distance(first, last) + 1);
    
    if (begins_before(mappable, *first)) {
        result.push_back(get_left_overhang(mappable, *first));
    }
    
    std::transform(first, std::prev(last), std::next(first), std::back_inserter(result),
                   [] (const auto& mappable, const auto& next_mappable) {
                       return get_intervening(mappable, next_mappable);
                   });
    
    if (ends_before(*std::prev(last), mappable)) {
        result.push_back(get_right_overhang(mappable, *std::prev(last)));
    }
    
    return result;
}

template <typename Container, typename Mappable>
std::vector<GenomicRegion> get_all_intervening(const Container& mappables, const Mappable& mappable)
{
    return get_all_intervening(std::cbegin(mappables), std::cend(mappables), mappable);
}

template <typename ForwardIterator>
auto segment_overlapped(ForwardIterator first, ForwardIterator last)
{
    std::vector<std::vector<typename ForwardIterator::value_type>> result {};
    
    if (first == last) return result;
    
    result.reserve(std::distance(first, last));
    
    auto it        = first;
    auto rightmost = first;
    
    while (first != last) {
        while (it != last && (overlaps(*it, *rightmost) || ends_equal(*it, *rightmost))) {
            if (ends_before(*rightmost, *it)) {
                rightmost = it;
            }
            ++it;
        }
        result.emplace_back(first, it);
        first = it;
        rightmost = it;
    }
    
    result.shrink_to_fit();
    
    return result;
}

template <typename Container>
auto segment_overlapped(const Container& mappables)
{
    return segment_overlapped(std::cbegin(mappables), std::cend(mappables));
}

template <typename ForwardIterator>
auto segment_by_begin(ForwardIterator first, ForwardIterator last)
{
    using MappableType = typename std::iterator_traits<ForwardIterator>::value_type;
    
    std::vector<std::vector<MappableType>> result {};
    
    if (first == last) return result;
    
    result.reserve(std::distance(first, last));
    
    while (first != last) {
        auto it = std::find_if_not(std::next(first), last, [first]
                                   (const auto& mappable) { return begins_equal(*first, mappable); });
        result.emplace_back(first, it);
        first = it;
    }
    
    result.shrink_to_fit();
    
    return result;
}

template <typename Container>
auto segment_by_begin(const Container& mappables)
{
    return segment_by_begin(std::cbegin(mappables), std::cend(mappables));
}

template <typename ForwardIterator>
auto segment_by_region(ForwardIterator first, ForwardIterator last)
{
    using MappableType = typename std::iterator_traits<ForwardIterator>::value_type;
    
    std::vector<std::vector<MappableType>> result {};
    
    if (first == last) return result;
    
    result.reserve(std::distance(first, last));
    
    while (first != last) {
        const auto& curr_region = get_region(*first);
        auto it = std::find_if_not(std::next(first), last, [&curr_region]
                                   (const auto& mappable) { return curr_region == get_region(mappable); });
        result.emplace_back(first, it);
        first = it;
    }
    
    result.shrink_to_fit();
    
    return result;
}

template <typename Container>
auto segment_by_region(const Container& mappables)
{
    return segment_equal(std::cbegin(mappables), std::cend(mappables));
}

template <typename Mappable>
std::vector<GenomicRegion> get_segment_regions(const std::vector<std::vector<Mappable>>& segments)
{
    std::vector<GenomicRegion> result {};
    result.reserve(segments.size());
    
    for (const auto& segment : segments) {
        result.push_back(get_encompassing_region(segment));
    }
    
    return result;
}

#endif
