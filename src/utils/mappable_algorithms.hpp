// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef mappable_algorithms_hpp
#define mappable_algorithms_hpp

#include <algorithm>
#include <numeric>
#include <cstddef>
#include <iterator>
#include <stdexcept>
#include <type_traits>

#include <boost/iterator/filter_iterator.hpp>
#include <boost/range/iterator_range_core.hpp>

#include "basics/contig_region.hpp"
#include "basics/genomic_region.hpp"
#include "concepts/mappable.hpp"
#include "concepts/mappable_range.hpp"
#include "type_tricks.hpp"

/**
 Mappable algorithms are STL like algorithms that work on ranges of Mappable objects.
 
 Most of the algorithms here require the input range to be sorted via Mappable::operator<
 (or meets the requirments of ForwardSorted).
 
 Some of the algorithms have lower time complexities when the input range also meets the
 requirement of BidirectionallySorted.
 */

namespace octopus {

// sum_region_sizes

template <typename InputIt>
auto sum_region_sizes(InputIt first, InputIt last)
{
    using MappableTp = typename std::iterator_traits<InputIt>::value_type;
    static_assert(is_region_or_mappable<MappableTp>, "Mappable required");
    using Position = typename RegionType<MappableTp>::Position;
    return std::accumulate(first, last, Position {0},
                           [] (const auto curr, const auto& mappable) {
                               return curr + region_size(mappable);
                           });
}

template <typename Container>
auto sum_region_sizes(const Container& mappables)
{
    return sum_region_sizes(std::cbegin(mappables), std::cend(mappables));
}

// leftmost_mappable

/**
 Returns the leftmost mappable element in the range [first, last).
 
 The range [first, last) is not required to be sorted.
 */
template <typename ForwardIt>
ForwardIt leftmost_mappable(ForwardIt first, ForwardIt last)
{
    using MappableTp = typename std::iterator_traits<ForwardIt>::value_type;
    static_assert(is_region_or_mappable<MappableTp>, "Mappable required");
    return std::min_element(first, last);
}

namespace detail {

template <typename C, typename = void>
struct HasMemberLeftmostMappable : std::false_type {};

template <typename C>
struct HasMemberLeftmostMappable<C, std::enable_if_t<
    std::is_same<decltype(std::declval<C>().leftmost_mappable()),
    typename C::const_iterator>::value>>
: std::true_type {};

template <typename Container>
auto leftmost_mappable(const Container& mappables, std::true_type)
{
    return mappables.leftmost_mappable();
}

template <typename Container>
auto leftmost_mappable(const Container& mappables, std::false_type)
{
    return leftmost_mappable(std::cbegin(mappables), std::cend(mappables));
}

} // namespace detail

template <typename Container>
auto leftmost_mappable(const Container& mappables)
{
    return detail::leftmost_mappable(mappables, detail::HasMemberLeftmostMappable<Container> {});
}

// rightmost_mappable

/**
 Returns the rightmost mappable element in the range [first, last).
 
 The range [first, last) is not required to be sorted.
 */
template <typename ForwardIt>
ForwardIt rightmost_mappable(ForwardIt first, ForwardIt last)
{
    using MappableTp = typename std::iterator_traits<ForwardIt>::value_type;
    static_assert(is_region_or_mappable<MappableTp>, "Mappable required");
    return std::max_element(first, last, [] (const auto& lhs, const auto& rhs) {
        return (ends_equal(lhs, rhs)) ? begins_before(lhs, rhs) : ends_before(lhs, rhs);
    });
}

namespace detail {

template <typename C, typename = void>
struct HasMemberRightmostMappable : std::false_type {};

template <typename C>
struct HasMemberRightmostMappable<C, std::enable_if_t<
    std::is_same<decltype(std::declval<C>().rightmost_mappable()),
    typename C::const_iterator>::value>>
: std::true_type {};

template <typename Container>
auto rightmost_mappable(const Container& mappables, std::true_type)
{
    return mappables.rightmost_mappable();
}

template <typename Container>
auto rightmost_mappable(const Container& mappables, std::false_type)
{
    return rightmost_mappable(std::cbegin(mappables), std::cend(mappables));
}

} // namespace detail

template <typename Container>
auto rightmost_mappable(const Container& mappables)
{
    return detail::rightmost_mappable(mappables, detail::HasMemberRightmostMappable<Container> {});
}

template <typename ForwardIt>
decltype(auto) leftmost_region(ForwardIt first, ForwardIt last)
{
    return mapped_region(*leftmost_mappable(first, last));
}

template <typename Container>
decltype(auto) leftmost_region(const Container& mappables)
{
    return mapped_region(*leftmost_mappable(mappables));
}

template <typename ForwardIt>
decltype(auto) rightmost_region(ForwardIt first, ForwardIt last)
{
    return mapped_region(*rightmost_mappable(first, last));
}

template <typename Container>
decltype(auto) rightmost_region(const Container& mappables)
{
    return mapped_region(*rightmost_mappable(mappables));
}

// largest_mappable

/**
 Returns the mappable element in the range [first, last) with the largest size.
 
 The range [first, last) is not required to be sorted.
 */
template <typename ForwardIt>
ForwardIt largest_mappable(ForwardIt first, ForwardIt last)
{
    using MappableTp = typename std::iterator_traits<ForwardIt>::value_type;
    static_assert(is_region_or_mappable<MappableTp>, "Mappable required");
    return std::max_element(first, last,
                            [] (const auto& lhs, const auto& rhs) {
                                return region_size(lhs) < region_size(rhs);
                            });
}

namespace detail {

template <typename C, typename = void>
struct HasMemberLargestMappable : std::false_type {};

template <typename C>
struct HasMemberLargestMappable<C, std::enable_if_t<
    std::is_same<decltype(std::declval<C>().largest_mappable()),
    typename C::const_iterator>::value>>
: std::true_type {};

template <typename Container>
auto largest_mappable(const Container& mappables, std::true_type)
{
    return mappables.largest_mappable();
}

template <typename Container>
auto largest_mappable(const Container& mappables, std::false_type)
{
    return largest_mappable(std::cbegin(mappables), std::cend(mappables));
}

} // namespace detail

template <typename Container>
auto largest_mappable(const Container& mappables)
{
    return detail::largest_mappable(mappables, detail::HasMemberLargestMappable<Container> {});
}

// smallest_mappable

/**
 Returns the mappable element in the range [first, last) with the smallest size.
 
 The range [first, last) is not required to be sorted.
 */
template <typename ForwardIt>
ForwardIt smallest_mappable(ForwardIt first, ForwardIt last)
{
    using MappableTp = typename std::iterator_traits<ForwardIt>::value_type;
    static_assert(is_region_or_mappable<MappableTp>, "Mappable required");
    return std::min_element(first, last,
                            [] (const auto& lhs, const auto& rhs) {
                                return region_size(lhs) < region_size(rhs);
                            });
}

namespace detail {

template <typename C, typename = void>
struct HasMemberSmallestMappable : std::false_type {};

template <typename C>
struct HasMemberSmallestMappable<C, std::enable_if_t<
std::is_same<decltype(std::declval<C>().smallest_mappable()),
    typename C::const_iterator>::value>>
: std::true_type {};

template <typename Container>
auto smallest_mappable(const Container& mappables, std::true_type)
{
    return mappables.smallest_mappable();
}

template <typename Container>
auto smallest_mappable(const Container& mappables, std::false_type)
{
    return smallest_mappable(std::cbegin(mappables), std::cend(mappables));
}

} // namespace detail

template <typename Container>
auto smallest_mappable(const Container& mappables)
{
    return detail::smallest_mappable(mappables, detail::HasMemberSmallestMappable<Container> {});
}

template <typename ForwardIt>
decltype(auto) largest_region(ForwardIt first, ForwardIt last)
{
    return mapped_region(*largest_mappable(first, last));
}

template <typename Container>
decltype(auto) largest_region(const Container& mappables)
{
    return mapped_region(*largest_mappable(mappables));
}

template <typename ForwardIt>
decltype(auto) smallest_region(ForwardIt first, ForwardIt last)
{
    return mapped_region(*smallest_mappable(first, last));
}

template <typename Container>
decltype(auto) smallest_region(const Container& mappables)
{
    return mapped_region(*smallest_mappable(mappables));
}

// is_bidirectionally_sorted

/**
 Returns true if the range of Mappable elements in the range [first, last) is meetis the
 requirments of BidirectionallySorted.
 */
template <typename ForwardIt>
bool is_bidirectionally_sorted(ForwardIt first, ForwardIt last)
{
    using MappableTp = typename std::iterator_traits<ForwardIt>::value_type;
    static_assert(is_region_or_mappable<MappableTp>, "Mappable required");
    return std::is_sorted(first, last,
                          [] (const auto& lhs, const auto& rhs) {
                              return lhs < rhs || ends_before(lhs, rhs);
                          });
}

template <typename Container>
bool is_bidirectionally_sorted(const Container& mappables)
{
    return is_bidirectionally_sorted(std::cbegin(mappables), std::cend(mappables));
}

// is_bidirectionally_sorted_until

/**
 Returns the an iterator to the first element in the range [first, last) that is not sorted
 according to BidirectionallySorted.
 */
template <typename ForwardIt>
ForwardIt is_bidirectionally_sorted_until(ForwardIt first, ForwardIt last)
{
    using MappableTp = typename std::iterator_traits<ForwardIt>::value_type;
    static_assert(is_region_or_mappable<MappableTp>, "Mappable required");
    return std::is_sorted_until(first, last,
                                [] (const auto& lhs, const auto& rhs) {
                                    return lhs < rhs || ends_before(lhs, rhs);
                                });
}

template <typename Container>
typename Container::const_iterator is_bidirectionally_sorted_until(const Container& mappables)
{
    return is_bidirectionally_sorted_until(std::cbegin(mappables), std::cend(mappables));
}

// extract_bidirectionally_sorted_ranges

/**
 Returns the minimum number of sub-ranges in each within the range [first, last) such that each
 sub-range is BidirectionallySorted.
 */
template <typename ForwardIt>
auto extract_bidirectionally_sorted_ranges(ForwardIt first, ForwardIt last)
{
    using MappableTp = typename std::iterator_traits<ForwardIt>::value_type;
    static_assert(is_region_or_mappable<MappableTp>, "Mappable required");
    std::vector<boost::iterator_range<ForwardIt>> result {};
    while (first != last) {
        auto it = is_bidirectionally_sorted_until(first, last);
        result.emplace_back(first, it);
        first = it;
    }
    result.shrink_to_fit();
    return result;
}

template <typename Container>
auto extract_bidirectionally_sorted_ranges(const Container& mappables)
{
    return extract_bidirectionally_sorted_ranges(std::cbegin(mappables), std::cend(mappables));
}

// find_first_after

/**
 Returns the first element in the range [first, last) that is_after mappable.
 
 Requires the range [first, last) is ForwardSorted.
 */
template <typename ForwardIt, typename MappableTp>
ForwardIt find_first_after(ForwardIt first, ForwardIt last, const MappableTp& mappable)
{
    using MappableTp2 = typename std::iterator_traits<ForwardIt>::value_type;
    static_assert(is_region_or_mappable<MappableTp> && is_region_or_mappable<MappableTp2>,
                  "Mappable required");
    auto it = std::lower_bound(first, last, next_mapped_position(mappable));
    return std::find_if_not(it, last, [&mappable] (const auto& m) { return overlaps(m, mappable); });
}

template <typename Container, typename MappableTp>
auto find_first_after(const Container& mappables, const MappableTp& mappable)
{
    return find_first_after(std::cbegin(mappables), std::cend(mappables), mappable);
}

// overlap_range

/**
 Returns an OverlapRange of the range [first, last).
 */
template <typename BidirIt, typename MappableTp>
OverlapRange<BidirIt> overlap_range(BidirIt first, BidirIt last, const MappableTp& mappable,
                                    ForwardSortedTag)
{
    using MappableTp2 = typename std::iterator_traits<BidirIt>::value_type;
    static_assert(is_region_or_mappable<MappableTp> && is_region_or_mappable<MappableTp2>,
                  "Mappable required");
    const auto it = find_first_after(first, last, mappable);
    // We must do a linear search for the first overlapped as the end positions may not be sorted.
    // Consider find the overlap range of M:
    //
    //    [--------A-------)
    //           [-----B------)
    //             [--M--)
    //
    // Here std::lower_bound (comparing mapped_begin) will return B and miss A.
    const auto it2 = std::find_if(first, it, [&mappable] (const auto& m) { return overlaps(m, mappable); });
    return make_overlap_range(it2, it, mappable);
}

template <typename BidirIt, typename MappableTp>
OverlapRange<BidirIt> overlap_range(BidirIt first, BidirIt last, const MappableTp& mappable,
                                    BidirectionallySortedTag)
{
    using MappableTp2 = typename std::iterator_traits<BidirIt>::value_type;
    static_assert(is_region_or_mappable<MappableTp> && is_region_or_mappable<MappableTp2>,
                  "Mappable required");
    auto overlapped = std::equal_range(first, last, mappable,
                                       [] (const auto& lhs, const auto& rhs) {
                                           return is_before(lhs, rhs);
                                       });
    // We need to try and push these boundaries out as the range does not fully capture
    // insertions
    overlapped.first = std::find_if_not(std::make_reverse_iterator(overlapped.first),
                                        std::make_reverse_iterator(first),
                                        [&mappable] (const auto& m) {
                                            return overlaps(m, mappable);
                                        }).base();
    overlapped.second = std::find_if_not(overlapped.second, last,
                                         [&mappable] (const auto& m) {
                                             return overlaps(m, mappable);
                                         });
    return make_overlap_range(overlapped.first, overlapped.second, mappable);
}

template <typename BidirIt, typename MappableTp>
auto overlap_range(BidirIt first, BidirIt last, const MappableTp& mappable)
{
    return overlap_range(first, last, mappable, ForwardSortedTag {});
}

// Faster version if the max mappable size in [first, last) is known
template <typename ForwardIt, typename MappableTp>
OverlapRange<ForwardIt>
overlap_range(ForwardIt first, ForwardIt last, const MappableTp& mappable,
              const typename RegionType<MappableTp>::Position max_mappable_size)
{
    using MappableTp2 = typename std::iterator_traits<ForwardIt>::value_type;
    static_assert(is_region_or_mappable<MappableTp> && is_region_or_mappable<MappableTp2>,
                  "Mappable required");
    const auto it1 = find_first_after(first, last, mappable);
    const auto leftmost = shift(mapped_region(mappable), -std::min(mapped_begin(mappable), max_mappable_size));
    auto it2 = std::lower_bound(first, it1, leftmost,
                                [] (const auto& lhs, const auto& rhs) {
                                    return begins_before(lhs, rhs);
                                });
    it2 = std::find_if(it2, it1, [&mappable] (const auto& m) { return overlaps(m, mappable); });
    return make_overlap_range(it2, it1, mappable);
}

namespace detail {

template <typename C, typename T, typename = void>
struct HasMemberOverlapRange : std::false_type {};

template <typename C, typename T>
struct HasMemberOverlapRange<C, T, std::enable_if_t<
std::is_same<decltype(std::declval<C>().overlap_range(std::declval<T>())),
    OverlapRange<typename C::const_iterator>>::value>>
: std::true_type {};

template <typename Container, typename MappableTp>
auto overlap_range(const Container& mappables, const MappableTp& mappable, std::true_type)
{
    return mappables.overlap_range(mappable);
}

template <typename Container, typename MappableTp>
auto overlap_range(const Container& mappables, const MappableTp& mappable, std::false_type)
{
    return overlap_range(std::cbegin(mappables), std::cend(mappables), mappable);
}

template <typename Container, typename MappableTp>
auto overlap_range(const Container& mappables, const MappableTp& mappable,
                   const typename RegionType<MappableTp>::Position max_mappable_size,
                   std::true_type)
{
    return mappables.overlap_range(mappable);
}

template <typename Container, typename MappableTp>
auto overlap_range(const Container& mappables, const MappableTp& mappable,
                   const typename RegionType<MappableTp>::Position max_mappable_size,
                   std::false_type)
{
    return overlap_range(std::cbegin(mappables), std::cend(mappables), mappable,
                           max_mappable_size);
}

} // namespace detail

template <typename Container, typename MappableTp>
auto overlap_range(const Container& mappables, const MappableTp& mappable)
{
    return detail::overlap_range(mappables, mappable,
                                 detail::HasMemberOverlapRange<Container, MappableTp> {});
}

template <typename Container, typename MappableTp>
OverlapRange<typename Container::const_iterator>
overlap_range(const Container& mappables, const MappableTp& mappable,
              const typename RegionType<MappableTp>::Position max_mappable_size)
{
    return detail::overlap_range(mappables, mappable, max_mappable_size,
                                 detail::HasMemberOverlapRange<Container, MappableTp> {});
}

template <typename Container, typename MappableTp>
OverlapRange<typename Container::const_iterator>
overlap_range(const Container& mappables, const MappableTp& mappable,
              BidirectionallySortedTag)
{
    return overlap_range(std::cbegin(mappables), std::cend(mappables), mappable,
                         BidirectionallySortedTag {});
}

// copy_overlapped

template <typename Container, typename MappableTp>
Container copy_overlapped(const Container& mappables, const MappableTp& mappable)
{
    const auto overlapped = overlap_range(mappables, mappable);
    return Container {std::cbegin(overlapped), std::cend(overlapped)};
}

template <typename Container, typename MappableTp>
Container copy_overlapped(const Container& mappables, const MappableTp& mappable,
                          BidirectionallySortedTag)
{
    const auto overlapped = overlap_range(mappables, mappable, BidirectionallySortedTag {});
    return Container {std::cbegin(overlapped), std::cend(overlapped)};
}

template <typename Container, typename MappableTp>
Container copy_overlapped(const Container& mappables, const MappableTp& mappable,
                          const typename RegionType<MappableTp>::Position max_mappable_size)
{
    const auto overlapped = overlap_range(mappables, mappable, max_mappable_size);
    return Container {std::cbegin(overlapped), std::cend(overlapped)};
}

// copy_nonoverlapped

template <typename Container, typename MappableType>
Container copy_nonoverlapped(const Container& mappables, const MappableType& mappable)
{
    using std::cbegin; using std::cend;
    const auto num_overlapped = count_overlapped(mappables, mappable);
    if (num_overlapped == 0) return mappables;
    Container result {};
    result.reserve(mappables.size() - num_overlapped);
    auto overlapped = overlap_range(mappables, mappable);
    auto base_begin = cbegin(overlapped).base();
    auto base_end   = cend(overlapped).base();
    result.insert(cbegin(mappables), base_begin);
    while (!overlapped.empty()) {
        overlapped.advance_begin(1);
        ++base_begin;
        if (overlapped.begin() != base_begin) {
            result.insert(base_begin, cbegin(overlapped).base());
            base_begin = cbegin(overlapped).base();
        }
    }
    result.insert(base_end, cend(mappables));
    return result;
}

// has_overlapped

/**
 Returns true if any of the mappable elements in the range [first, last) overlap with mappable.
 
 Requires the range [first, last) is ForwardSorted.
 */
template <typename BidirIt, typename MappableTp>
bool has_overlapped(BidirIt first, BidirIt last, const MappableTp& mappable,
                    ForwardSortedTag)
{
    using MappableTp2 = typename std::iterator_traits<BidirIt>::value_type;
    static_assert(is_region_or_mappable<MappableTp> && is_region_or_mappable<MappableTp2>,
                  "Mappable required");
    if (first == last) return false;
    const auto it = find_first_after(first, last, mappable);
    const auto it2 = std::find_if_not(it, last, [&mappable] (const auto& m) { return overlaps(m, mappable); });
    if (it != it2) return true;
    // searches in reverse order on the assumption regions closer to the boundry with
    // mappable are more likely to overlap with mappable.
    return std::any_of(std::make_reverse_iterator(it), std::make_reverse_iterator(first),
                       [&mappable] (const auto& m) {
                           return overlaps(mappable, m);
                       });
}

template <typename BidirIt, typename MappableTp>
bool has_overlapped(BidirIt first, BidirIt last, const MappableTp& mappable,
                    BidirectionallySortedTag)
{
    using MappableTp2 = typename std::iterator_traits<BidirIt>::value_type;
    static_assert(is_region_or_mappable<MappableTp> && is_region_or_mappable<MappableTp2>,
                  "Mappable required");
    return std::binary_search(first, last, mappable,
                              [] (const auto& lhs, const auto& rhs) { return is_before(lhs, rhs); });
}

template <typename BidirIt, typename MappableTp>
bool has_overlapped(BidirIt first, BidirIt last, const MappableTp& mappable,
                    const typename RegionType<MappableTp>::Position max_mappable_size)
{
    const auto overlapped = overlap_range(first, last, mappable, max_mappable_size);
    return !overlapped.empty();
}

template <typename BidirIt, typename MappableTp>
bool has_overlapped(BidirIt first, BidirIt last, const MappableTp& mappable)
{
    return has_overlapped(first, last, mappable, ForwardSortedTag {});
}

namespace detail {

template <typename C, typename T, typename = void>
struct HasMemberHasOverlapped : std::false_type {};

template <typename C, typename T>
struct HasMemberHasOverlapped<C, T, std::enable_if_t<
std::is_same<decltype(std::declval<C>().has_overlapped(std::declval<T>())), bool>::value>>
: std::true_type {};

template <typename Container, typename MappableTp>
auto has_overlapped(const Container& mappables, const MappableTp& mappable, std::true_type)
{
    return mappables.has_overlapped(mappable);
}

template <typename Container, typename MappableTp>
auto has_overlapped(const Container& mappables, const MappableTp& mappable, std::false_type)
{
    return has_overlapped(std::cbegin(mappables), std::cend(mappables), mappable);
}

template <typename Container, typename MappableTp>
auto has_overlapped(const Container& mappables, const MappableTp& mappable,
                   const typename RegionType<MappableTp>::Position max_mappable_size,
                   std::true_type)
{
    return mappables.has_overlapped(mappable);
}

template <typename Container, typename MappableTp>
auto has_overlapped(const Container& mappables, const MappableTp& mappable,
                    const typename RegionType<MappableTp>::Position max_mappable_size,
                    std::false_type)
{
    return has_overlapped(std::cbegin(mappables), std::cend(mappables), mappable,
                          max_mappable_size);
}

} // namespace detail

template <typename Container, typename MappableTp>
auto has_overlapped(const Container& mappables, const MappableTp& mappable)
{
    return detail::has_overlapped(mappables, mappable,
                                 detail::HasMemberHasOverlapped<Container, MappableTp> {});
}

template <typename Container, typename MappableTp>
auto
has_overlapped(const Container& mappables, const MappableTp& mappable,
               BidirectionallySortedTag)
{
    return has_overlapped(std::cbegin(mappables), std::cend(mappables), mappable,
                          BidirectionallySortedTag {});
}

template <typename Container, typename MappableTp>
auto
has_overlapped(const Container& mappables, const MappableTp& mappable,
               const typename RegionType<MappableTp>::Position max_mappable_size)
{
    return detail::has_overlapped(mappables, mappable, max_mappable_size,
                                  detail::HasMemberHasOverlapped<Container, MappableTp> {});
}

// count_overlapped

/**
 Returns the number of mappable elements in the range [first, last) that overlap with mappable.
 
 Requires the range [first, last) is ForwardSorted.
 */
template <typename ForwardIt, typename MappableTp>
std::size_t count_overlapped(ForwardIt first, ForwardIt last, const MappableTp& mappable,
                             ForwardSortedTag)
{
    using MappableTp2 = typename std::iterator_traits<ForwardIt>::value_type;
    static_assert(is_region_or_mappable<MappableTp> && is_region_or_mappable<MappableTp2>,
                  "Mappable required");
    const auto overlapped = overlap_range(first, last, mappable);
    return size(overlapped, ForwardSortedTag {});
}

template <typename ForwardIt, typename MappableTp>
std::size_t count_overlapped(ForwardIt first, ForwardIt last, const MappableTp& mappable,
                             BidirectionallySortedTag)
{
    using MappableTp2 = typename std::iterator_traits<ForwardIt>::value_type;
    static_assert(is_region_or_mappable<MappableTp> && is_region_or_mappable<MappableTp2>,
                  "Mappable required");
    const auto overlapped = overlap_range(first, last, mappable, BidirectionallySortedTag {});
    return size(overlapped, BidirectionallySortedTag {});
}

/**
 Returns the number of mappable elements in the range [first, last) that overlap with mappable.
 
 Requires the range [first, last) is ForwardSorted.
 */
template <typename ForwardIt, typename MappableTp>
std::size_t count_overlapped(ForwardIt first, ForwardIt last, const MappableTp& mappable,
                             GenomicRegion::Position max_mappable_size)
{
    using MappableTp2 = typename std::iterator_traits<ForwardIt>::value_type;
    static_assert(is_region_or_mappable<MappableTp> && is_region_or_mappable<MappableTp2>,
                  "Mappable required");
    const auto overlapped = overlap_range(first, last, mappable, max_mappable_size);
    return std::distance(std::cbegin(overlapped), std::cend(overlapped));
}

template <typename ForwardIt, typename MappableTp>
std::size_t count_overlapped(ForwardIt first, ForwardIt last, const MappableTp& mappable)
{
    return count_overlapped(first, last, mappable, ForwardSortedTag {});
}

namespace detail {

template <typename C, typename T, typename = void>
struct HasMemberCountOverlapped : std::false_type {};

template <typename C, typename T>
struct HasMemberCountOverlapped<C, T, std::enable_if_t<
std::is_same<decltype(std::declval<C>().count_overlapped(std::declval<T>())), std::size_t>::value>>
: std::true_type {};

template <typename Container, typename MappableTp>
auto count_overlapped(const Container& mappables, const MappableTp& mappable, std::true_type)
{
    return mappables.count_overlapped(mappable);
}

template <typename Container, typename MappableTp>
auto count_overlapped(const Container& mappables, const MappableTp& mappable, std::false_type)
{
    return count_overlapped(std::cbegin(mappables), std::cend(mappables), mappable);
}

template <typename Container, typename MappableTp>
auto count_overlapped(const Container& mappables, const MappableTp& mappable,
                    const typename RegionType<MappableTp>::Position max_mappable_size,
                    std::true_type)
{
    return mappables.count_overlapped(mappable);
}

template <typename Container, typename MappableTp>
auto count_overlapped(const Container& mappables, const MappableTp& mappable,
                    const typename RegionType<MappableTp>::Position max_mappable_size,
                    std::false_type)
{
    return count_overlapped(std::cbegin(mappables), std::cend(mappables), mappable,
                            max_mappable_size);
}

} // namespace detail

template <typename Container, typename MappableTp>
auto count_overlapped(const Container& mappables, const MappableTp& mappable)
{
    return detail::count_overlapped(mappables, mappable,
                                  detail::HasMemberCountOverlapped<Container, MappableTp> {});
}

template <typename Container, typename MappableTp>
auto count_overlapped(const Container& mappables, const MappableTp& mappable,
                      const typename RegionType<MappableTp>::Position max_mappable_size)
{
    return detail::count_overlapped(mappables, mappable, max_mappable_size,
                                  detail::HasMemberCountOverlapped<Container, MappableTp> {});
}

// has_exact_overlap

/**
 Returns true if the range [first, last) has an element with the region of mappable
 
 Requires the range [first, last) is ForwardSorted.
 */
template <typename ForwardIt, typename MappableTp>
bool has_exact_overlap(ForwardIt first, ForwardIt last, const MappableTp& mappable,
                       ForwardSortedTag)
{
    using MappableTp2 = typename std::iterator_traits<ForwardIt>::value_type;
    static_assert(is_region_or_mappable<MappableTp> && is_region_or_mappable<MappableTp2>,
                  "Mappable required");
    const auto overlapped = overlap_range(first, last, mappable);
    return std::find_if(std::cbegin(overlapped), std::cend(overlapped),
                        [&mappable] (const auto& e) {
                            return is_same_region(mappable, e);
                        }) != std::cend(overlapped);
}

template <typename ForwardIt, typename MappableTp>
bool has_exact_overlap(ForwardIt first, ForwardIt last, const MappableTp& mappable,
                       BidirectionallySortedTag)
{
    using MappableTp2 = typename std::iterator_traits<ForwardIt>::value_type;
    static_assert(is_region_or_mappable<MappableTp> && is_region_or_mappable<MappableTp2>,
                  "Mappable required");
    return std::binary_search(first, last, mappable);
}

template <typename ForwardIt, typename MappableTp>
bool has_exact_overlap(ForwardIt first, ForwardIt last, const MappableTp& mappable)
{
    return has_exact_overlap(first, last, mappable, ForwardSortedTag {});
}

template <typename Container, typename MappableTp>
bool has_exact_overlap(const Container& container, const MappableTp& mappable,
                       ForwardSortedTag)
{
    return has_exact_overlap(std::cbegin(container), std::cend(container), mappable,
                             ForwardSortedTag {});
}

template <typename Container, typename MappableTp>
bool has_exact_overlap(const Container& container, const MappableTp& mappable,
                       BidirectionallySortedTag)
{
    return has_exact_overlap(std::cbegin(container), std::cend(container), mappable,
                             BidirectionallySortedTag {});
}

template <typename Container, typename MappableTp>
bool has_exact_overlap(const Container& container, const MappableTp& mappable)
{
    return has_exact_overlap(container, mappable, ForwardSortedTag {});
}

// contained_range

/**
 Returns an ContainedRange of the range [first, last).
 
 Requires the range [first, last) is ForwardSorted.
 */
template <typename BidirIt, typename MappableTp>
ContainedRange<BidirIt>
contained_range(BidirIt first, BidirIt last, const MappableTp& mappable)
{
    using MappableTp2 = typename std::iterator_traits<BidirIt>::value_type;
    static_assert(is_region_or_mappable<MappableTp> && is_region_or_mappable<MappableTp2>,
                  "Mappable required");
    const auto it = std::lower_bound(first, last, mappable,
                               [] (const auto& lhs, const auto& rhs) {
                                   return begins_before(lhs, rhs);
                               });
    const auto it2 = find_first_after(it, last, mappable);
    if (it == it2) return make_contained_range(it, it2, mappable);
    auto rit = std::find_if(std::make_reverse_iterator(it2), std::make_reverse_iterator(std::next(it)),
                            [&mappable] (const auto& m) { return contains(mappable, m); });
    return make_contained_range(it, rit.base(), mappable);
}

namespace detail {

template <typename C, typename T, typename = void>
struct HasMemberContainedRange : std::false_type {};

template <typename C, typename T>
struct HasMemberContainedRange<C, T, std::enable_if_t<
    std::is_same<decltype(std::declval<C>().contained_range(std::declval<T>())),
        ContainedRange<typename C::const_iterator>>::value>>
: std::true_type {};

template <typename Container, typename MappableTp>
auto contained_range(const Container& mappables, const MappableTp& mappable, std::true_type)
{
    return mappables.contained_range(mappable);
}

template <typename Container, typename MappableTp>
auto contained_range(const Container& mappables, const MappableTp& mappable, std::false_type)
{
    return contained_range(std::cbegin(mappables), std::cend(mappables), mappable);
}

} // namespace detail

template <typename Container, typename MappableTp>
auto contained_range(const Container& mappables, const MappableTp& mappable)
{
    return detail::contained_range(mappables, mappable,
                                    detail::HasMemberContainedRange<Container, MappableTp> {});
}

// has_contained

/**
 Returns true if any of the mappable elements in the range [first, last) are contained within mappable.
 
 Requires the range [first, last) is ForwardSorted.
 */
template <typename ForwardIt, typename MappableTp>
bool has_contained(ForwardIt first, ForwardIt last, const MappableTp& mappable)
{
    using MappableTp2 = typename std::iterator_traits<ForwardIt>::value_type;
    static_assert(is_region_or_mappable<MappableTp> && is_region_or_mappable<MappableTp2>,
                  "Mappable required");
    auto it = std::lower_bound(first, last, mappable,
                               [] (const auto& lhs, const auto& rhs) {
                                   return begins_before(lhs, rhs);
                               });
    return (it != last) && mapped_end(*it) <= mapped_end(mappable);
}

namespace detail {

template <typename C, typename T, typename = void>
struct HasMemberHasContained : std::false_type {};

template <typename C, typename T>
struct HasMemberHasContained<C, T, std::enable_if_t<
    std::is_same<decltype(std::declval<C>().has_contained(std::declval<T>())), bool>::value>>
: std::true_type {};

template <typename Container, typename MappableTp>
auto has_contained(const Container& mappables, const MappableTp& mappable, std::true_type)
{
    return mappables.has_contained(mappable);
}

template <typename Container, typename MappableTp>
auto has_contained(const Container& mappables, const MappableTp& mappable, std::false_type)
{
    return has_contained(std::cbegin(mappables), std::cend(mappables), mappable);
}

} // namespace detail

template <typename Container, typename MappableTp>
auto has_contained(const Container& mappables, const MappableTp& mappable)
{
    return detail::has_contained(mappables, mappable,
                                   detail::HasMemberHasContained<Container, MappableTp> {});
}

// count_contained

/**
 Returns the number of mappable elements in the range [first, last) that are contained within mappable.
 
 Requires the range [first, last) is ForwardSorted.
 */
template <typename BidirIt, typename MappableTp>
std::size_t count_contained(BidirIt first, BidirIt last, const MappableTp& mappable)
{
    using MappableTp2 = typename std::iterator_traits<BidirIt>::value_type;
    static_assert(is_region_or_mappable<MappableTp> && is_region_or_mappable<MappableTp2>,
                  "Mappable required");
    const auto contained = contained_range(first, last, mappable);
    return std::distance(contained.begin(), contained.end());
}

namespace detail {

template <typename C, typename T, typename = void>
struct HasMemberCountContained : std::false_type {};

template <typename C, typename T>
struct HasMemberCountContained<C, T, std::enable_if_t<
std::is_same<decltype(std::declval<C>().count_contained(std::declval<T>())), std::size_t>::value>>
: std::true_type {};

template <typename Container, typename MappableTp>
auto count_contained(const Container& mappables, const MappableTp& mappable, std::true_type)
{
    return mappables.count_contained(mappable);
}

template <typename Container, typename MappableTp>
auto count_contained(const Container& mappables, const MappableTp& mappable, std::false_type)
{
    return count_contained(std::cbegin(mappables), std::cend(mappables), mappable);
}

} // namespace detail

template <typename Container, typename MappableTp>
auto count_contained(const Container& mappables, const MappableTp& mappable)
{
    return detail::count_contained(mappables, mappable,
                                 detail::HasMemberCountContained<Container, MappableTp> {});
}

// copy_contained

template <typename Container, typename MappableType>
Container copy_contained(const Container& mappables, const MappableType& mappable)
{
    const auto contained = contained_range(mappables, mappable);
    return Container {std::begin(contained), std::end(contained)};
}

// copy_contained

template <typename Container, typename MappableType>
Container copy_noncontained(const Container& mappables, const MappableType& mappable)
{
    using std::cbegin; using std::cend;
    const auto num_overlapped = count_overlapped(mappables, mappable);
    if (num_overlapped == 0) return mappables;
    Container result {};
    result.reserve(mappables.size() - num_overlapped);
    auto contained  = contained_range(mappables, mappable);
    auto base_begin = cbegin(contained).base();
    auto base_end   = cend(contained).base();
    result.insert(cbegin(mappables), base_begin);
    while (!contained.empty()) {
        contained.advance_begin(1);
        std::advance(base_begin);
        if (contained.begin() != base_begin) {
            result.insert(base_begin, cbegin(contained).base());
            base_begin = cbegin(contained).base();
        }
    }
    result.insert(base_end, cend(mappables));
    return result;
}

// count_spanning

template <typename Container, typename MappableTp>
std::size_t count_spanning(const Container& mappables, const MappableTp& mappable)
{
    using MappableTp2 = typename Container::value_type;
    static_assert(is_region_or_mappable<MappableTp> && is_region_or_mappable<MappableTp2>,
                  "Mappable required");
    const auto overlapped = overlap_range(mappables, mappable);
    return std::count_if(std::cbegin(overlapped), std::cend(overlapped),
                         [&mappable] (const auto& m) {
                             return contains(mapped_region(m), mappable);
                         });
}

// count_shared

/**
 Returns the number of Mappable elements in the range [first, last) that both lhs and rhs overlap.
 
 Requires the range [first, last) is ForwardSorted.
 */
template <typename BidirIt, typename MappableTp1, typename MappableTp2, typename OrderTag>
std::size_t count_shared(BidirIt first, BidirIt last,
                         const MappableTp1& lhs, const MappableTp2& rhs,
                         OrderTag)
{
    using MappableTp3 = typename std::iterator_traits<BidirIt>::value_type;
    static_assert(is_region_or_mappable<MappableTp1> && is_region_or_mappable<MappableTp2>
                  && is_region_or_mappable<MappableTp3>,
                  "Mappable required");
    const auto lhs_overlapped = overlap_range(first, last, lhs, OrderTag {});
    const auto rhs_overlapped = overlap_range(first, last, rhs, OrderTag {});
    return (size(lhs_overlapped) <= size(rhs_overlapped)) ?
            std::count_if(lhs_overlapped.begin(), lhs_overlapped.end(),
                          [&rhs] (const auto& region) { return overlaps(region, rhs); }) :
            std::count_if(rhs_overlapped.begin(), rhs_overlapped.end(),
                          [&lhs] (const auto& region) { return overlaps(region, lhs); });
}

template <typename Container, typename MappableTp1, typename MappableTp2, typename OrderTag>
auto count_shared(const Container& container, const MappableTp1& lhs, const MappableTp2& rhs,
                  OrderTag)
{
    return count_shared(std::cbegin(container), std::cend(container), lhs, rhs, OrderTag {});
}

/**
 Returns if any of the Mappable elements in the range [first, last) overlaps both lhs and rhs.
 
 Requires the range [first, last) is ForwardSorted.
 */
template <typename BidirIt, typename MappableTp1, typename MappableTp2, typename OrderTag>
bool has_shared(BidirIt first, BidirIt last,
                const MappableTp1& lhs, const MappableTp2& rhs,
                OrderTag)
{
    using MappableTp3 = typename std::iterator_traits<BidirIt>::value_type;
    static_assert(is_region_or_mappable<MappableTp1> && is_region_or_mappable<MappableTp2>
                  && is_region_or_mappable<MappableTp3>,
                  "Mappable required");
    const auto lhs_overlapped = overlap_range(first, last, lhs, OrderTag {});
    const auto rhs_overlapped = overlap_range(first, last, rhs, OrderTag {});
    return (size(lhs_overlapped) <= size(rhs_overlapped)) ?
            std::any_of(lhs_overlapped.begin(), lhs_overlapped.end(),
                        [&rhs] (const auto& region) { return overlaps(region, rhs); }) :
            std::any_of(rhs_overlapped.begin(), rhs_overlapped.end(),
                        [&lhs] (const auto& region) { return overlaps(region, lhs); });
}

template <typename Container, typename MappableTp1, typename MappableTp2, typename OrderTag>
bool has_shared(const Container& container, const MappableTp1& lhs, const MappableTp2& rhs,
                OrderTag)
{
    return has_shared(std::cbegin(container), std::cend(container), lhs, rhs, OrderTag {});
}

/**
 Returns the first Mappable element in the range [first2, last2) such that the element shares a region
 in the range [first1, last1) with mappable.
 
 Requires [first1, last1) and [first2, last2) are sorted w.r.t GenomicRegion::operator<
 */
template <typename BidirIt1, typename BidirIt2, typename MappableTp, typename OrderTag>
BidirIt2 find_first_shared(BidirIt1 first1, BidirIt1 last1,
                           BidirIt2 first2, BidirIt2 last2,
                           const MappableTp& mappable,
                           OrderTag)
{
    using MappableTp2 = typename std::iterator_traits<BidirIt1>::value_type;
    using MappableTp3 = typename std::iterator_traits<BidirIt2>::value_type;
    static_assert(is_region_or_mappable<MappableTp> && is_region_or_mappable<MappableTp2>
                  && is_region_or_mappable<MappableTp3>,
                  "Mappable required");
    return std::find_if(first2, last2,
                        [first1, last1, &mappable] (const auto& m) {
                            return has_shared(first1, last1, m, mappable, OrderTag {});
                        });
}

/**
 Counts the number of Mappable elements in the range [first2 + 1, last2) that have shared
 elements in the range [first1, last1) with first2
 
 Requires [first1, last1) and [first2, last2) to be sorted w.r.t GenomicRegion::operator<
 */
template <typename BidirIt1, typename BidirIt2, typename OrderTag>
std::size_t count_if_shared_with_first(BidirIt1 first1, BidirIt1 last1,
                                       BidirIt2 first2, BidirIt2 last2,
                                       OrderTag)
{
    using MappableTp1 = typename std::iterator_traits<BidirIt1>::value_type;
    using MappableTp2 = typename std::iterator_traits<BidirIt2>::value_type;
    static_assert(is_region_or_mappable<MappableTp1> && is_region_or_mappable<MappableTp2>,
                  "Mappable required");
    if (first2 == last2) return 0;
    const auto overlapped = overlap_range(first1, last1, *first2, OrderTag {});
    if (empty(overlapped)) return 0;
    return size(overlap_range(std::next(first2), last2, *std::prev(overlapped.end()), OrderTag {}));
}

template <typename Container, typename ForwardIt>
std::size_t count_if_shared_with_first(const Container& mappables,
                                       ForwardIt first, ForwardIt last)
{
    if (first == last) return 0;
    const auto overlapped = overlap_range(mappables, *first);
    if (empty(overlapped)) return 0;
    return count_overlapped(std::next(first), last, overlapped.back());
}

// adjacent_overlap_find

/**
 Returns the first element of the first adjacent pair of elements in the range [first, last) that
 overlap, or last if non are found.
 */
template <typename ForwardIt>
ForwardIt adjacent_overlap_find(ForwardIt first, const ForwardIt last)
{
    using MappableTp = typename std::iterator_traits<ForwardIt>::value_type;
    static_assert(is_region_or_mappable<MappableTp>, "Mappable required");
    if (first == last) return last;
    auto prev = first++;
    while (first != last) {
        if (overlaps(*prev, *first++)) return prev;
        ++prev;
    }
    return last;
}

template <typename Container>
auto adjacent_overlap_find(const Container& mappables)
{
    return adjacent_overlap_find(std::cbegin(mappables), std::cend(mappables));
}

// extract_regions

namespace detail {

template <typename InputIt>
auto extract_regions(InputIt first, InputIt last, std::input_iterator_tag)
{
    using MappableTp = typename std::iterator_traits<InputIt>::value_type;
    static_assert(is_region_or_mappable<MappableTp>, "Mappable required");
    std::vector<RegionType<MappableTp>> result {};
    std::transform(first, last, std::back_inserter(result),
                   [] (const auto& mappable) { return mapped_region(mappable); });
    result.shrink_to_fit();
    return result;
}

template <typename ForwardIt>
auto extract_regions(ForwardIt first, ForwardIt last, std::forward_iterator_tag)
{
    using MappableTp = typename std::iterator_traits<ForwardIt>::value_type;
    static_assert(is_region_or_mappable<MappableTp>, "Mappable required");
    std::vector<RegionType<MappableTp>> result {};
    result.reserve(std::distance(first, last));
    std::transform(first, last, std::back_inserter(result),
                   [] (const auto& mappable) { return mapped_region(mappable); });
    return result;
}

} // namespace detail

/**
 Returns a vector of mapped_regions in the range [first, last).
 */
template <typename InputIt>
auto extract_regions(InputIt first, InputIt last)
{
    using Category = typename std::iterator_traits<InputIt>::iterator_category;
    return detail::extract_regions(first, last, Category {});
}

template <typename Container>
auto extract_regions(const Container& mappables)
{
    return extract_regions(std::cbegin(mappables), std::cend(mappables));
}

// decompose

namespace detail {

template <typename MappableTp>
auto decompose(const MappableTp& mappable, ContigRegion)
{
    std::vector<ContigRegion> result {};
    const auto num_elements = size(mappable);
    if (num_elements == 0) return result;
    result.reserve(num_elements);
    ContigRegion::Position n {0};
    std::generate_n(std::back_inserter(result), num_elements, [&] () {
        const auto begin = mapped_begin(mappable) + n;
        ++n;
        return ContigRegion {begin, begin + 1};
    });
    return result;
}

template <typename MappableTp>
auto decompose(const MappableTp& mappable, GenomicRegion)
{
    std::vector<GenomicRegion> result {};
    const auto num_elements = region_size(mappable);
    if (num_elements == 0) return result;
    result.reserve(num_elements);
    GenomicRegion::Position n {0};
    const auto& contig = contig_name(mappable);
    std::generate_n(std::back_inserter(result), num_elements, [&] () {
        const auto begin = mapped_begin(mappable) + n;
        ++n;
        return GenomicRegion {contig, begin, begin + 1};
    });
    return result;
}

} // namespace detail

/**
 Returns a vector of RegionType<MappableTp>'s, each of size 1, that cover the region defined
 by mappable.
 */
template <typename MappableTp>
auto decompose(const MappableTp& mappable)
{
    static_assert(is_region_or_mappable<MappableTp>, "Mappable required");
    return detail::decompose(mappable, RegionType<MappableTp> {});
}

/**
 Returns the maximal vector of RegionType<MappableTp>'s, each of size n, that do not
 span past mapped_end(mappable).
 */
template <typename MappableTp>
auto decompose(const MappableTp& mappable, const GenomicRegion::Position n)
{
    static_assert(is_region_or_mappable<MappableTp>, "Mappable required");
    std::vector<GenomicRegion> result {};
    if (n == 0) return result;
    const auto num_elements = region_size(mappable) / n;
    if (num_elements == 0) return result;
    result.reserve(num_elements);
    const auto& contig = contig_name(mappable);
    auto curr = mapped_begin(mappable);
    std::generate_n(std::back_inserter(result), num_elements, [&contig, &curr, n] () {
        auto tmp = curr;
        curr += n;
        return GenomicRegion {contig, tmp, tmp + n};
    });
    return result;
}

// encompassing_region

/**
 Returns the region encompassed by the elements in the range [first, last).
 
 Requires the range [first, last) is ForwardSorted.
 */
template <typename ForwardIt, typename = enable_if_iterator<ForwardIt>>
auto encompassing_region(ForwardIt first, ForwardIt last)
{
    using MappableTp = typename std::iterator_traits<ForwardIt>::value_type;
    static_assert(is_region_or_mappable<MappableTp>, "Mappable required");
    if (first == last) {
        throw std::logic_error {"get_encompassing given empty range"};
    }
    return encompassing_region(*first, *rightmost_mappable(first, last));
}

template <typename Container>
auto encompassing_region(const Container& mappables)
{
    return encompassing_region(leftmost_region(mappables), rightmost_region(mappables));
}

// extract_*_regions

namespace detail {

template <typename ForwardIt, typename Compare>
auto extract_overlapping_regions(ForwardIt first, const ForwardIt last, Compare cmp)
{
    using MappableTp = typename std::iterator_traits<ForwardIt>::value_type;
    static_assert(is_region_or_mappable<MappableTp>, "Mappable required");
    std::vector<RegionType<MappableTp>> result {};
    if (first == last) return result;
    result.reserve(std::distance(first, last));
    auto first_overlapped = first;
    auto rightmost        = first;
    for (; first != last; ++first) {
        if (cmp(mapped_begin(*first), mapped_end(*rightmost))) {
            if (result.empty() || !ends_equal(result.back(), *rightmost)) {
                result.push_back(closed_region(*first_overlapped, *rightmost));
            }
            rightmost        = first;
            first_overlapped = first;
        } else if (ends_before(*rightmost, *first)) {
            rightmost = first;
        }
    }
    result.push_back(closed_region(*first_overlapped, *rightmost));
    result.shrink_to_fit();
    return result;
}

} // namespace detail

// extract_covered_regions

/**
 Returns the minimal range of non-overlapping regions such that each element in the range [first, last)
 is contained within a single region.
 
 Requires [first_mappable, last_mappable) is sorted w.r.t GenomicRegion::operator<
 */
template <typename ForwardIt>
auto extract_covered_regions(ForwardIt first, ForwardIt last)
{
    return detail::extract_overlapping_regions(first, last,
                                               [] (const auto& a, const auto b) {
                                                   return a > b;
                                               });
}

template <typename Container>
auto extract_covered_regions(const Container& mappables)
{
    return extract_covered_regions(std::cbegin(mappables), std::cend(mappables));
}

// extract_mutually_exclusive_regions

/**
 Returns the maximal range of non-overlapping regions such that each element in the range [first, last)
 is contained within a single region.
 
 Requires [first_mappable, last_mappable) is sorted w.r.t GenomicRegion::operator<
 */

template <typename ForwardIt>
auto extract_mutually_exclusive_regions(ForwardIt first, const ForwardIt last)
{
    return detail::extract_overlapping_regions(first, last,
                                               [] (const auto& a, const auto b) {
                                                   return a >= b;
                                               });
}

template <typename Container>
auto extract_mutually_exclusive_regions(const Container& mappables)
{
    return extract_mutually_exclusive_regions(std::cbegin(mappables), std::cend(mappables));
}

// extract_intervening_regions

/**
 Returns all intervening regions between non-overlapping mappables in the range [first, last).
 
 Requires the range [first, last) is ForwardSorted.
 */
template <typename ForwardIt>
auto extract_intervening_regions(ForwardIt first, ForwardIt last)
{
    using MappableTp = typename std::iterator_traits<ForwardIt>::value_type;
    static_assert(is_region_or_mappable<MappableTp>, "Mappable required");
    std::vector<RegionType<MappableTp>> result {};
    if (first == last) return result;
    result.reserve(std::distance(first, last) - 1);
    std::transform(first, std::prev(last), std::next(first), std::back_inserter(result),
                   [] (const auto& mappable, const auto& next_mappable) {
                       return *intervening_region(mappable, next_mappable);
                   });
    return result;
}

template <typename Container>
auto extract_intervening_regions(const Container& mappables)
{
    return extract_intervening_regions(std::cbegin(mappables), std::cend(mappables));
}

/**
 Returns all intervening regions between non-overlapping mappables in the range
 [first, last), and also any flanking regions of mappable if the range [first, last) is 
 contained within mappable.
 
 Requires the range [first, last) is ForwardSorted.
 */
template <typename ForwardIt, typename MappableTp>
auto extract_intervening_regions(ForwardIt first, ForwardIt last, const MappableTp& mappable)
{
    using MappableTp2 = typename std::iterator_traits<ForwardIt>::value_type;
    static_assert(is_region_or_mappable<MappableTp> && is_region_or_mappable<MappableTp2>,
                  "Mappable required");
    std::vector<RegionType<MappableTp>> result {};
    if (first == last) return result;
    result.reserve(std::distance(first, last) + 1);
    if (begins_before(mappable, *first)) {
        result.push_back(left_overhang_region(mappable, *first));
    }
    std::transform(first, std::prev(last), std::next(first), std::back_inserter(result),
                   [] (const auto& mappable, const auto& next_mappable) {
                       return *intervening_region(mappable, next_mappable);
                   });
    if (ends_before(*std::prev(last), mappable)) {
        result.push_back(right_overhang_region(mappable, *std::prev(last)));
    }
    return result;
}

template <typename Container, typename Mappable>
auto extract_intervening_regions(const Container& mappables, const Mappable& mappable)
{
    return extract_intervening_regions(std::cbegin(mappables), std::cend(mappables), mappable);
}

// segment_*

template <typename ForwardIt>
auto segment_overlapped_copy(ForwardIt first, ForwardIt last)
{
    using MappableTp = typename std::iterator_traits<ForwardIt>::value_type;
    static_assert(is_region_or_mappable<MappableTp>, "Mappable required");
    std::vector<std::vector<MappableTp>> result {};
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
        first     = it;
        rightmost = it;
    }
    result.shrink_to_fit();
    return result;
}

template <typename Container>
auto segment_overlapped_copy(const Container& mappables)
{
    return segment_overlapped_copy(std::cbegin(mappables), std::cend(mappables));
}

template <typename Container>
auto segment_by_overlapped_move(Container& mappables)
{
    return segment_overlapped_copy(std::make_move_iterator(std::begin(mappables)),
                                   std::make_move_iterator(std::end(mappables)));
}

template <typename ForwardIt>
auto segment_by_begin_copy(ForwardIt first, ForwardIt last)
{
    using MappableTp = typename std::iterator_traits<ForwardIt>::value_type;
    static_assert(is_region_or_mappable<MappableTp>, "Mappable required");
    std::vector<std::vector<MappableTp>> result {};
    if (first == last) return result;
    result.reserve(std::distance(first, last));
    while (first != last) {
        auto it = std::find_if_not(std::next(first), last, [first]
                                   (const auto& mappable) {
                                       return begins_equal(*first, mappable);
                                   });
        result.emplace_back(first, it);
        first = it;
    }
    result.shrink_to_fit();
    return result;
}

template <typename Container>
auto segment_by_begin_copy(const Container& mappables)
{
    return segment_by_begin_copy(std::cbegin(mappables), std::cend(mappables));
}

template <typename Container>
auto segment_by_begin_move(Container& mappables)
{
    return segment_by_begin_move(std::make_move_iterator(std::begin(mappables)),
                                 std::make_move_iterator(std::end(mappables)));
}

template <typename ForwardIt>
auto segment_by_end_copy(ForwardIt first, ForwardIt last)
{
    using MappableTp = typename std::iterator_traits<ForwardIt>::value_type;
    static_assert(is_region_or_mappable<MappableTp>, "Mappable required");
    std::vector<std::vector<MappableTp>> result {};
    if (first == last) return result;
    result.reserve(std::distance(first, last));
    while (first != last) {
        auto it = std::find_if_not(std::next(first), last, [first]
                                   (const auto& mappable) {
                                       return ends_equal(*first, mappable);
                                   });
        result.emplace_back(first, it);
        first = it;
    }
    result.shrink_to_fit();
    return result;
}

template <typename Container>
auto segment_by_end_copy(const Container& mappables)
{
    return segment_by_end_copy(std::cbegin(mappables), std::cend(mappables));
}

template <typename Container>
auto segment_by_end_move(Container& mappables)
{
    return segment_by_end_copy(std::make_move_iterator(std::begin(mappables)),
                               std::make_move_iterator(std::end(mappables)));
}

template <typename ForwardIt>
auto segment_by_region_copy(ForwardIt first, ForwardIt last)
{
    using MappableTp = typename std::iterator_traits<ForwardIt>::value_type;
    static_assert(is_region_or_mappable<MappableTp>, "Mappable required");
    std::vector<std::vector<MappableTp>> result {};
    if (first == last) return result;
    result.reserve(std::distance(first, last));
    while (first != last) {
        const auto& curr_region = mapped_region(*first);
        auto it = std::find_if_not(std::next(first), last, [&curr_region]
                                   (const auto& mappable) {
                                       return curr_region == mapped_region(mappable);
                                   });
        result.emplace_back(first, it);
        first = it;
    }
    result.shrink_to_fit();
    return result;
}

template <typename Container>
auto segment_by_region_copy(const Container& mappables)
{
    return segment_equal(std::cbegin(mappables), std::cend(mappables));
}

template <typename MappableTp>
auto all_segment_regions(const std::vector<std::vector<MappableTp>>& segments)
{
    static_assert(is_region_or_mappable<MappableTp>, "Mappable required");
    std::vector<RegionType<MappableTp>> result {};
    result.reserve(segments.size());
    for (const auto& segment : segments) {
        result.push_back(encompassing_region(segment));
    }
    return result;
}

// calculate_positional_coverage

/**
 Returns the number of elements that overlap each position within region.
 */
template <typename ForwardIt, typename RegionTp,
          typename = EnableIfRegionOrMappable<typename std::iterator_traits<ForwardIt>::value_type>>
auto calculate_positional_coverage(ForwardIt first, ForwardIt last, const RegionTp& region)
{
    using MappableTp = typename std::iterator_traits<ForwardIt>::value_type;
    static_assert(std::is_same<RegionType<MappableTp>, RegionTp>::value,
                  "RegionType mismatch");
    using std::next; using std::min;
    const auto num_positions = region_size(region);
    std::vector<unsigned> result(num_positions, 0);
    const auto result_begin_itr = std::begin(result);
    const auto first_position = mapped_begin(region);
    std::for_each(first, last, [=] (const auto& mappable) {
        const auto it1 = next(result_begin_itr,
                              (mapped_begin(mappable) <= first_position)
                              ? 0 : mapped_begin(mappable) - first_position);
        const auto it2 = next(result_begin_itr, min(mapped_end(mappable) - first_position, num_positions));
        std::transform(it1, it2, it1, [] (const auto count) { return count + 1; });
    });
    return result;
}

template <typename ForwardIt,
          typename = EnableIfRegionOrMappable<typename std::iterator_traits<ForwardIt>::value_type>>
auto calculate_positional_coverage(ForwardIt first, ForwardIt last)
{
    return calculate_positional_coverage(first, last, encompassing_region(first, last));
}

template <typename Container,
          typename = EnableIfMappable<typename Container::value_type>>
auto calculate_positional_coverage(const Container& mappables)
{
    return calculate_positional_coverage(std::cbegin(mappables), std::cend(mappables));
}

template <typename Container, typename RegionTp,
          typename = EnableIfRegionOrMappable<typename Container::value_type>>
auto calculate_positional_coverage(const Container& mappables, const RegionTp& region)
{
    using MappableTp = typename Container::value_type;
    static_assert(std::is_same<RegionType<MappableTp>, RegionTp>::value,
                  "RegionType mismatch");
    const auto overlapped = overlap_range(mappables, region);
    return calculate_positional_coverage(std::cbegin(overlapped), std::cend(overlapped), region);
}

template <typename Container, typename RegionTp,
          typename = enable_if_not_map<Container>,
          typename = EnableIfRegionOrMappable<typename Container::value_type>>
bool has_coverage(const Container& mappables, const RegionTp& region)
{
    using MappableTp = typename Container::value_type;
    static_assert(std::is_same<RegionType<MappableTp>, RegionTp>::value,
                  "RegionType mismatch");
    if (mappables.empty() || is_empty(region)) return false;
    const auto overlapped = overlap_range(mappables, region);
    return std::any_of(std::cbegin(overlapped), std::cend(overlapped),
                       [] (const auto& mappable) {
                           return !is_empty_region(mappable);
                       });
}

template <typename Container,
          typename = enable_if_not_map<Container>,
          typename = EnableIfRegionOrMappable<typename Container::value_type>>
bool has_coverage(const Container& mappables)
{
    return std::any_of(std::cbegin(mappables), std::cend(mappables),
                       [] (const auto& mappable) {
                           return !is_empty_region(mappable);
                       });
}

template <typename Container, typename RegionTp,
          typename = enable_if_not_map<Container>,
          typename = EnableIfRegionOrMappable<typename Container::value_type>>
unsigned min_coverage(const Container& mappables, const RegionTp& region)
{
    using MappableTp = typename Container::value_type;
    static_assert(std::is_same<RegionType<MappableTp>, RegionTp>::value,
                  "RegionType mismatch");
    if (mappables.empty() || is_empty(region)) return 0;
    const auto positional_coverage = calculate_positional_coverage(mappables, region);
    return *std::min_element(std::cbegin(positional_coverage), std::cend(positional_coverage));
}

template <typename Container,
          typename = enable_if_not_map<Container>,
          typename = EnableIfRegionOrMappable<typename Container::value_type>>
unsigned min_coverage(const Container& mappables)
{
    if (mappables.empty()) return 0;
    const auto positional_coverage = calculate_positional_coverage(mappables);
    return *std::min_element(std::cbegin(positional_coverage), std::cend(positional_coverage));
}

template <typename Container, typename RegionTp,
          typename = enable_if_not_map<Container>,
          typename = EnableIfRegionOrMappable<typename Container::value_type>>
unsigned max_coverage(const Container& mappables, const RegionTp& region)
{
    using MappableTp = typename Container::value_type;
    static_assert(std::is_same<RegionType<MappableTp>, RegionTp>::value,
                  "RegionType mismatch");
    if (mappables.empty() || is_empty(region)) return 0;
    const auto positional_coverage = calculate_positional_coverage(mappables, region);
    return *std::max_element(std::cbegin(positional_coverage), std::cend(positional_coverage));
}

template <typename Container,
          typename = enable_if_not_map<Container>,
          typename = EnableIfRegionOrMappable<typename Container::value_type>>
unsigned max_coverage(const Container& mappables)
{
    if (mappables.empty()) return 0;
    const auto positional_coverage = calculate_positional_coverage(mappables);
    return *std::max_element(std::cbegin(positional_coverage), std::cend(positional_coverage));
}

} // namespace octopus

#endif
