//
//  mappable_algorithms.hpp
//  Octopus
//
//  Created by Daniel Cooke on 10/04/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_mappable_algorithms_hpp
#define Octopus_mappable_algorithms_hpp

#include <algorithm>
#include <numeric>
#include <cstddef>
#include <iterator>
#include <stdexcept>
#include <type_traits>

#include <boost/iterator/filter_iterator.hpp>
#include <boost/range/iterator_range_core.hpp>

#include "contig_region.hpp"
#include "genomic_region.hpp"
#include "mappable.hpp"
#include "mappable_ranges.hpp"

#include <iostream> // DEBUG

template <typename ForwardIt>
auto sum_sizes(ForwardIt first, ForwardIt last)
{
    return std::accumulate(first, last, GenomicRegion::SizeType {},
                           [] (const auto curr, const auto& mappable) { return curr + size(mappable); });
}

template <typename Container>
auto sum_sizes(const Container& mappables)
{
    return sum_sizes(std::cbegin(mappables), std::cend(mappables));
}

/**
 Returns the leftmost mappable element in the range [first, last).
 
 The range [first, last) is not required to be sorted.
 */
template <typename ForwardIt>
ForwardIt leftmost_mappable(ForwardIt first, ForwardIt last)
{
    return std::min_element(first, last);
}

template <typename Container>
auto leftmost_mappable(const Container& mappables)
{
    return leftmost_mappable(std::cbegin(mappables), std::cend(mappables));
}

/**
 Returns the rightmost mappable element in the range [first, last).
 
 The range [first, last) is not required to be sorted.
 */
template <typename ForwardIt>
inline
ForwardIt rightmost_mappable(ForwardIt first, ForwardIt last)
{
    return std::max_element(first, last, [] (const auto& lhs, const auto& rhs) {
        return (ends_equal(lhs, rhs)) ? begins_before(lhs, rhs) : ends_before(lhs, rhs);
    });
}

template <typename Container>
auto rightmost_mappable(const Container& mappables)
{
    return rightmost_mappable(std::cbegin(mappables), std::cend(mappables));
}

/**
 Returns the mappable element in the range [first, last) with the largest size.
 
 The range [first, last) is not required to be sorted.
 */
template <typename ForwardIt>
inline
ForwardIt largest_element(ForwardIt first, ForwardIt last)
{
    return std::max_element(first, last,
                            [] (const auto& lhs, const auto& rhs) { return size(lhs) < size(rhs); });
}

template <typename Container>
auto largest_element(const Container& mappables)
{
    return largest_element(std::cbegin(mappables), std::cend(mappables));
}

/**
 Returns the mappable element in the range [first, last) with the smallest size.
 
 The range [first, last) is not required to be sorted.
 */
template <typename ForwardIt>
inline
ForwardIt smallest_element(ForwardIt first, ForwardIt last)
{
    return std::min_element(first, last,
                            [] (const auto& lhs, const auto& rhs) { return size(lhs) < size(rhs); });
}

template <typename Container>
auto smallest_element(const Container& mappables)
{
    return smallest_element(std::cbegin(mappables), std::cend(mappables));
}

/**
 Returns the first Mappable element in the range [first, last) that is_after mappable
 
 Requires [first, last) is sorted w.r.t GenomicRegion::operator<
 */
template <typename ForwardIt, typename MappableTp>
inline
ForwardIt find_first_after(ForwardIt first, ForwardIt last, const MappableTp& mappable)
{
    return std::lower_bound(first, last, next_position(mappable));
}

template <typename Container, typename MappableTp>
auto find_first_after(const Container& mappables, const MappableTp& mappable)
{
    return find_first_after(std::cbegin(mappables), std::cend(mappables), mappable);
}

/**
 Returns true if the range of Mappable elements in the range [first, last) is sorted
 w.r.t GenomicRegion::operator<, and satisfies the condition 'if lhs < rhs then end(lhs) < end(rhs)
 */
template <typename ForwardIt>
bool is_bidirectionally_sorted(ForwardIt first, ForwardIt last)
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
template <typename ForwardIt>
ForwardIt is_bidirectionally_sorted_until(ForwardIt first, ForwardIt last)
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
template <typename ForwardIt>
std::vector<boost::iterator_range<ForwardIt>>
bidirectionally_sorted_ranges(ForwardIt first, ForwardIt last)
{
    std::vector<boost::iterator_range<ForwardIt>> result {};
    
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
template <typename ForwardIt, typename MappableTp>
OverlapRange<ForwardIt>
overlap_range(ForwardIt first, ForwardIt last, const MappableTp& mappable,
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

template <typename Container, typename MappableTp>
OverlapRange<typename Container::const_iterator>
overlap_range(const Container& mappables, const MappableTp& mappable,
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
template <typename ForwardIt, typename MappableTp>
OverlapRange<ForwardIt>
overlap_range(ForwardIt first, ForwardIt last, const MappableTp& mappable,
              GenomicRegion::SizeType max_mappable_size)
{
    using MappableTp2 = typename std::iterator_traits<ForwardIt>::value_type;
    
    auto it = find_first_after(first, last, mappable);
    
    auto it2 = std::lower_bound(first, it, shift(get_region(mappable),
                                                 -std::min(get_begin(mappable), max_mappable_size)),
                                [] (const auto& lhs, const auto& rhs) { return begins_before(lhs, rhs); });
    
    it2 = std::find_if(it2, it, [&mappable] (const auto& m) { return overlaps(m, mappable); });
    
    return make_overlap_range(it2, it, mappable);
}

template <typename Container, typename MappableTp>
OverlapRange<typename Container::const_iterator>
overlap_range(const Container& mappables, const MappableTp& mappable,
              GenomicRegion::SizeType max_mappable_size)
{
    return overlap_range(std::cbegin(mappables), std::cend(mappables), mappable, max_mappable_size);
}

template <typename Container, typename MappableTp>
Container
copy_overlapped(const Container& mappables, const MappableTp& mappable,
                MappableRangeOrder order = MappableRangeOrder::ForwardSorted)
{
    const auto overlapped = overlap_range(mappables, mappable, order);
    return Container {std::cbegin(overlapped), std::cend(overlapped)};
}

template <typename Container, typename MappableTp>
Container
copy_overlapped(const Container& mappables, const MappableTp& mappable,
                GenomicRegion::SizeType max_mappable_size)
{
    const auto overlapped = overlap_range(mappables, mappable, max_mappable_size);
    return Container {std::cbegin(overlapped), std::cend(overlapped)};
}

/**
 Returns true if any of the mappable elements in the range [first, last) overlap with mappable.
 
 Requires [first, last) is sorted w.r.t GenomicRegion::operator<
 */
template <typename BidirectionalIterator, typename MappableTp>
bool has_overlapped(BidirectionalIterator first, BidirectionalIterator last,
                    const MappableTp& mappable,
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

template <typename Container, typename MappableTp>
bool has_overlapped(const Container& container, const MappableTp& mappable,
                    MappableRangeOrder order = MappableRangeOrder::ForwardSorted)
{
    return has_overlapped(std::cbegin(container), std::cend(container), mappable, order);
}

/**
 Returns the number of mappable elements in the range [first, last) that overlap with mappable.
 
 Requires [first, last) is sorted w.r.t GenomicRegion::operator<
 */
template <typename ForwardIt, typename MappableTp>
size_t count_overlapped(ForwardIt first, ForwardIt last, const MappableTp& mappable,
                        MappableRangeOrder order = MappableRangeOrder::ForwardSorted)
{
    const auto overlapped = overlap_range(first, last, mappable, order);
    return std::distance(overlapped.begin(), overlapped.end());
}

template <typename Container, typename MappableTp>
size_t count_overlapped(const Container& container, const MappableTp& mappable,
                        MappableRangeOrder order = MappableRangeOrder::ForwardSorted)
{
    return count_overlapped(std::cbegin(container), std::cend(container), mappable, order);
}

/**
 Returns the number of mappable elements in the range [first, last) that overlap with mappable.
 
 Requires [first, last) is sorted w.r.t GenomicRegion::operator<
 */
template <typename ForwardIt, typename MappableTp>
size_t count_overlapped(ForwardIt first, ForwardIt last, const MappableTp& mappable,
                        GenomicRegion::SizeType max_mappable_size)
{
    const auto overlapped = overlap_range(first, last, mappable, max_mappable_size);
    return std::distance(overlapped.begin(), overlapped.end());
}

template <typename Container, typename MappableTp>
size_t count_overlapped(const Container& container, const MappableTp& mappable,
                        GenomicRegion::SizeType max_mappable_size)
{
    return count_overlapped(std::cbegin(container), std::cend(container), mappable, max_mappable_size);
}

/**
 Returns true if the range [first, last) contains the exact region given my mappable
 
 Requires [first, last) is sorted w.r.t GenomicRegion::operator<
 */
template <typename ForwardIt, typename MappableTp>
bool has_exact_overlap(ForwardIt first, ForwardIt last, const MappableTp& mappable,
                       MappableRangeOrder order = MappableRangeOrder::ForwardSorted)
{
    if (order == MappableRangeOrder::BidirectionallySorted) {
        return std::binary_search(first, last, mappable);
    }
    
    const auto overlapped = overlap_range(first, last, mappable);
    
    return std::find_if(std::cbegin(overlapped), std::cend(overlapped),
                        [&mappable] (const auto& e) { return is_same_region(mappable, e); }) != std::cend(overlapped);
}

template <typename Container, typename MappableTp>
bool has_exact_overlap(const Container& container, const MappableTp& mappable,
                       MappableRangeOrder order = MappableRangeOrder::ForwardSorted)
{
    return has_exact_overlap(std::cbegin(container), std::cend(container), mappable, order);
}

/**
 Returns the sub-range of Mappable elements in the range [first, last) such that each element
 in the sub-range is contained within mappable.
 
 Requires [first, last) is sorted w.r.t GenomicRegion::operator<
 */
template <typename BidirectionalIterator, typename MappableTp>
ContainedRange<BidirectionalIterator>
contained_range(BidirectionalIterator first, BidirectionalIterator last, const MappableTp& mappable)
{
    auto it = std::lower_bound(first, last, mappable,
                               [] (const auto& lhs, const auto& rhs) { return begins_before(lhs, rhs); });
    
    auto it2 = find_first_after(it, last, mappable);
    
    if (it == it2) return make_contained_range(it, it2, mappable);
    
    auto rit = std::find_if(std::make_reverse_iterator(it2), std::make_reverse_iterator(std::next(it)),
                            [&mappable] (const auto& m) { return contains(mappable, m); });
    
    return make_contained_range(it, rit.base(), mappable);
}

template <typename Container, typename MappableTp>
ContainedRange<typename Container::const_iterator>
contained_range(const Container& mappables, const MappableTp& mappable)
{
    return contained_range(std::cbegin(mappables), std::cend(mappables), mappable);
}

/**
 Returns true if any of the mappable elements in the range [first, last) are contained within mappable.
 
 Requires [first, last) is sorted w.r.t GenomicRegion::operator<
 */
template <typename BidirectionalIterator, typename MappableTp>
bool has_contained(BidirectionalIterator first, BidirectionalIterator last, const MappableTp& mappable)
{
    auto it = std::lower_bound(first, last, mappable,
                               [] (const auto& lhs, const auto& rhs) { return begins_before(lhs, rhs); });
    
    return (it != last) && get_end(*it) <= get_end(mappable);
}

template <typename Container, typename MappableTp>
bool has_contained(const Container& container, const MappableTp& mappable)
{
    return has_contained(std::cbegin(container), std::cend(container), mappable);
}

/**
 Returns the number of mappable elements in the range [first, last) that are contained within mappable.
 
 Requires [first, last) is sorted w.r.t GenomicRegion::operator<
 */
template <typename BidirectionalIterator, typename MappableTp>
size_t count_contained(BidirectionalIterator first, BidirectionalIterator last, const MappableTp& mappable)
{
    const auto contained = contained_range(first, last, mappable);
    return std::distance(contained.begin(), contained.end());
}

template <typename Container, typename MappableTp>
size_t count_contained(const Container& container, const MappableTp& mappable)
{
    return count_contained(std::cbegin(container), std::cend(container), mappable);
}

/**
 Returns the number of Mappable elements in the range [first, last) that both lhs and rhs overlap.
 
 Requires [first, last) is sorted w.r.t GenomicRegion::operator<
 */
template <typename ForwardIt, typename MappableTp1, typename MappableTp2>
size_t count_shared(ForwardIt first, ForwardIt last,
                    const MappableTp1& lhs, const MappableTp2& rhs,
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

template <typename Container, typename MappableTp1, typename MappableTp2>
size_t count_shared(const Container& container,
                       const MappableTp1& lhs, const MappableTp2& rhs,
                       MappableRangeOrder order = MappableRangeOrder::ForwardSorted)
{
    return count_shared(std::cbegin(container), std::cend(container), lhs, rhs, order);
}

template <typename ForwardIt>
bool has_exact(ForwardIt first, ForwardIt last, const typename ForwardIt::value_type& mappable)
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
template <typename ForwardIt, typename MappableTp1, typename MappableTp2>
bool has_shared(ForwardIt first, ForwardIt last,
                const MappableTp1& lhs, const MappableTp2& rhs,
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
template <typename ForwardIt1, typename ForwardIt2, typename MappableTp>
ForwardIt2 find_first_shared(ForwardIt1 first1, ForwardIt1 last1,
                                   ForwardIt2 first2, ForwardIt2 last2,
                                   const MappableTp& mappable,
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
template <typename ForwardIt1, typename ForwardIt2>
size_t count_if_shared_with_first(ForwardIt1 first1, ForwardIt1 last1,
                                  ForwardIt2 first2, ForwardIt2 last2,
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
    std::vector<RegionType<typename Container::value_type>> result {};
    result.reserve(mappables.size());
    
    for (const auto& mappable : mappables) {
        result.emplace_back(get_region(mappable));
    }
    
    return result;
}

namespace detail
{
    template <typename MappableTp>
    auto decompose(const MappableTp& mappable, ContigRegion)
    {
        std::vector<ContigRegion> result {};
        
        const auto num_elements = size(mappable);
        
        if (num_elements == 0) return result;
        
        result.reserve(num_elements);
        
        ContigRegion::SizeType n {0};
        
        std::generate_n(std::back_inserter(result), num_elements, [&] () {
            return ContigRegion {get_begin(mappable) + n, get_begin(mappable) + ++n};
        });
        
        return result;
    }
    
    template <typename MappableTp>
    auto decompose(const MappableTp& mappable, GenomicRegion)
    {
        std::vector<GenomicRegion> result {};
        
        const auto num_elements = size(mappable);
        
        if (num_elements == 0) return result;
        
        result.reserve(num_elements);
        
        GenomicRegion::SizeType n {0};
        
        const auto& contig = get_contig_name(mappable);
        
        std::generate_n(std::back_inserter(result), num_elements, [&] () {
            return GenomicRegion {contig, get_begin(mappable) + n, get_begin(mappable) + ++n};
        });
        
        return result;
    }
} // namespace detail

/**
 Splits mappable into an ordered vector of GenomicRegions of size 1
 */
template <typename MappableTp>
auto decompose(const MappableTp& mappable)
{
    return detail::decompose(mappable, RegionType<MappableTp> {});
}

template <typename MappableTp>
auto decompose(const MappableTp& mappable, GenomicRegion::SizeType n)
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
template <typename ForwardIt>
auto get_encompassing_region(ForwardIt first, ForwardIt last)
{
    if (first == last) {
        throw std::runtime_error {"get_encompassing_region given empty range"};
    }
    
    return get_encompassing(*first, *rightmost_mappable(first, last));
}

template <typename Container>
auto get_encompassing_region(const Container& mappables)
{
    return get_encompassing_region(std::cbegin(mappables), std::cend(mappables));
}

/**
 Returns the minimal range of non-overlapping GenomicRegion's such that each element in the range [first, last)
 is contained within a single region.
 
 Requires [first_mappable, last_mappable) is sorted w.r.t GenomicRegion::operator<
 */
template <typename ForwardIt>
auto get_covered_regions(ForwardIt first, ForwardIt last)
{
    using ResultType = std::vector<RegionType<typename std::iterator_traits<ForwardIt>::value_type>>;
    
    ResultType result {};
    
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
std::vector<GenomicRegion>
get_covered_regions(const Container& mappables)
{
    return get_covered_regions(std::cbegin(mappables), std::cend(mappables));
}

/**
 Returns all intervening GenomicRegion's between non-overlapping mappables in the range 
 [first, last)
 
 Requires [first, last) is sorted w.r.t GenomicRegion::operator<
 */
template <typename ForwardIt>
auto get_all_intervening(ForwardIt first, ForwardIt last)
{
    using ResultType = std::vector<RegionType<typename std::iterator_traits<ForwardIt>::value_type>>;
    
    if (first == last) return ResultType {};
    
    ResultType result {};
    result.reserve(std::distance(first, last) - 1);
    
    std::transform(first, std::prev(last), std::next(first), std::back_inserter(result),
                   [] (const auto& mappable, const auto& next_mappable) {
                       return get_intervening(mappable, next_mappable);
                   });
    
    return result;
}

template <typename Container>
std::vector<GenomicRegion>
get_all_intervening(const Container& mappables)
{
    return get_all_intervening(std::cbegin(mappables), std::cend(mappables));
}

/**
 Returns all intervening GenomicRegion's between non-overlapping mappables in the range
 [first, last), and also any flanking regions of mappable if the range [first, last) is 
 contained within mappable.
 
 Requires [first, last) is sorted w.r.t GenomicRegion::operator<
 */
template <typename ForwardIt, typename Mappable>
auto get_all_intervening(ForwardIt first, ForwardIt last, const Mappable& mappable)
{
    using ResultType = std::vector<RegionType<Mappable>>;
    
    static_assert(std::is_same<RegionType<typename std::iterator_traits<ForwardIt>::value_type>,
                  RegionType<Mappable>>::value,
                  "RegionType of input range must match RegionType of mappable");
    
    if (first == last) return ResultType {};
    
    ResultType result {};
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
auto get_all_intervening(const Container& mappables, const Mappable& mappable)
{
    return get_all_intervening(std::cbegin(mappables), std::cend(mappables), mappable);
}

template <typename ForwardIt>
auto segment_overlapped(ForwardIt first, ForwardIt last)
{
    std::vector<std::vector<typename ForwardIt::value_type>> result {};
    
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

template <typename ForwardIt>
auto segment_by_begin(ForwardIt first, ForwardIt last)
{
    using MappableTp = typename std::iterator_traits<ForwardIt>::value_type;
    
    std::vector<std::vector<MappableTp>> result {};
    
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

template <typename ForwardIt>
auto segment_by_region(ForwardIt first, ForwardIt last)
{
    using MappableTp = typename std::iterator_traits<ForwardIt>::value_type;
    
    std::vector<std::vector<MappableTp>> result {};
    
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
auto get_segment_regions(const std::vector<std::vector<Mappable>>& segments)
{
    std::vector<RegionType<Mappable>> result {};
    result.reserve(segments.size());
    
    for (const auto& segment : segments) {
        result.push_back(get_encompassing_region(segment));
    }
    
    return result;
}

//template <typename InputIt>
//bool has_coverage(InputIt first, InputIt last, const GenomicRegion& region)
//{
//    const auto overlapped = overlap_range(first, last);
//    return std::any_of(std::cbegin(overlapped), std::cend(overlapped),
//                       [] (const auto& mappable) { return !empty(mappable); });
//}
//
//template <typename InputIt>
//bool has_coverage(InputIt first, InputIt last)
//{
//    return std::any_of(first, last, [] (const auto& mappable) { return !empty(mappable); });
//}
//
//template <typename Container>
//bool has_coverage(const Container& mappables)
//{
//    return has_coverage(std::cbegin(mappables), std::cend(mappables));
//}
//
//template <typename Container>
//bool has_coverage(const Container& mappables, const GenomicRegion& region)
//{
//    return has_coverage(std::cbegin(mappables), std::cend(mappables), region);
//}

template <typename InputIt, typename RegionTp>
auto positional_coverage(InputIt first, InputIt last, const RegionTp& region)
{
    static_assert(std::is_same<RegionType<typename std::iterator_traits<InputIt>::value_type>,
                  RegionTp>::value,
                  "RegionType of input range must match RegionTp");
    
    using std::for_each; using std::next; using std::min;
    
    const auto num_positions = size(region);
    
    std::vector<unsigned> result(num_positions, 0);
    
    const auto result_begin_itr = std::begin(result);
    
    const auto first_position = get_begin(region);
    
    for_each(first, last, [result_begin_itr, first_position, num_positions] (const auto& mappable) {
        const auto it1 = next(result_begin_itr, (get_begin(mappable) <= first_position) ? 0 : get_begin(mappable) - first_position);
        const auto it2 = next(result_begin_itr, min(get_end(mappable) - first_position, num_positions));
        std::transform(it1, it2, it1, [] (const auto count) { return count + 1; });
    });
    
    return result;
}

template <typename Container>
auto positional_coverage(const Container& mappables)
{
    return positional_coverage(std::cbegin(mappables), std::cend(mappables),
                               get_encompassing_region(mappables));
}

template <typename Container, typename RegionTp>
auto positional_coverage(const Container& mappables, const RegionTp& region)
{
    static_assert(std::is_same<RegionType<typename Container::value_type>, RegionTp>::value,
                  "RegionType of input container must match RegionTp");
    
    const auto overlapped = overlap_range(mappables, region);
    return positional_coverage(std::cbegin(overlapped), std::cend(overlapped), region);
}

#endif
