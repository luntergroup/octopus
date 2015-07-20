//
//  read_utils.h
//  Octopus
//
//  Created by Daniel Cooke on 09/04/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__read_utils__
#define __Octopus__read_utils__

#include <vector>
#include <unordered_map>
#include <iterator>  // std::begin, std::end, std::make_move_iterator
#include <utility>   // std::move
#include <algorithm> // std::min_element, std::max_element, std::transform, std::for_each

#include "aligned_read.h"
#include "read_filter.h"
#include "context_iterators.h"
#include "mappable_algorithms.h"
#include "maths.h"
#include "mappable_set.h"
#include "mappable_map.h"

#include <iostream> // TEST

using Octopus::ContextInserter;
using Octopus::ContextBackInserter;

namespace detail {
    template <typename Container>
    void reserve_if_enabled(Container& container, typename Container::size_type n) {}
    
    template <typename T>
    void reserve_if_enabled(std::vector<T>& container, typename std::vector<T>::size_type n) { container.reserve(n); }
    
    template <typename Container>
    void shrink_to_fit_if_enabled(Container& container) {}
    
    template <typename T>
    void shrink_to_fit_if_enabled(std::vector<T>& container) { container.shrink_to_fit(); }
}

template <typename ReadFilter>
std::pair<MappableSet<AlignedRead>, MappableSet<AlignedRead>>
filter_reads(const MappableSet<AlignedRead>& reads, ReadFilter& filter)
{
    MappableSet<AlignedRead> good_reads {}, bad_reads {};
    good_reads.reserve(reads.size());
    bad_reads.reserve(reads.size());
    
    filter.filter_reads(std::cbegin(reads), std::cend(reads), ContextInserter(good_reads), ContextInserter(bad_reads));
    
    good_reads.shrink_to_fit();
    bad_reads.shrink_to_fit();
    
    return {good_reads, bad_reads};
}

template <typename ReadFilter>
std::pair<MappableSet<AlignedRead>, MappableSet<AlignedRead>>
filter_reads(MappableSet<AlignedRead>&& reads, ReadFilter& filter)
{
    MappableSet<AlignedRead> good_reads {}, bad_reads {};
    good_reads.reserve(reads.size());
    bad_reads.reserve(reads.size());
    
    filter.filter_reads(std::make_move_iterator(std::begin(reads)), std::make_move_iterator(std::end(reads)),
                        ContextInserter(good_reads), ContextInserter(bad_reads));
    
    good_reads.shrink_to_fit();
    bad_reads.shrink_to_fit();
    
    return {good_reads, bad_reads};
}

template <typename KeyType, typename ReadFilter>
std::pair<MappableMap<KeyType, AlignedRead>, MappableMap<KeyType, AlignedRead>>
filter_reads(const MappableMap<KeyType, AlignedRead>& reads, ReadFilter& filter)
{
    MappableMap<KeyType, AlignedRead> good_reads {}, bad_reads {};
    good_reads.reserve(reads.size());
    bad_reads.reserve(reads.size());
    
    for (const auto& sample_reads : reads) {
        auto sample_filtered_reads = filter_reads(sample_reads.second, filter);
        good_reads.emplace(sample_reads.first, std::move(sample_filtered_reads.first));
        bad_reads.emplace(sample_reads.first, std::move(sample_filtered_reads.second));
    }
    
    return {good_reads, bad_reads};
}

template <typename KeyType, typename ReadFilter>
std::pair<MappableMap<KeyType, AlignedRead>, MappableMap<KeyType, AlignedRead>>
filter_reads(MappableMap<KeyType, AlignedRead>&& reads, ReadFilter& filter)
{
    MappableMap<KeyType, AlignedRead> good_reads {}, bad_reads {};
    good_reads.reserve(reads.size());
    bad_reads.reserve(reads.size());
    
    for (const auto& sample_reads : reads) {
        auto sample_filtered_reads = filter_reads(std::move(sample_reads.second), filter);
        good_reads.emplace(sample_reads.first, std::move(sample_filtered_reads.first));
        bad_reads.emplace(sample_reads.first, std::move(sample_filtered_reads.second));
    }
    
    return {good_reads, bad_reads};
}

template <typename T, typename Container>
using ReadMap = std::unordered_map<T, Container>;

template <typename T, typename Container, typename ReadFilter>
std::pair<ReadMap<T, Container>, ReadMap<T, Container>>
filter_reads(ReadMap<T, Container>&& reads, ReadFilter& read_filter)
{
    ReadMap<T, Container> good_read_map {}, bad_read_map {};
    good_read_map.reserve(reads.size());
    bad_read_map.reserve(reads.size());
    
    for (auto& sample_reads : reads) {
        Container good_reads {}, bad_reads {};
        
        detail::reserve_if_enabled(good_reads, sample_reads.second.size());
        detail::reserve_if_enabled(bad_reads, sample_reads.second.size() / 10); // arbitrarily chosen
        
        read_filter.filter_reads(std::make_move_iterator(std::begin(sample_reads.second)),
                                 std::make_move_iterator(std::end(sample_reads.second)),
                                 ContextBackInserter(good_reads), ContextBackInserter(bad_reads));
        
        sample_reads.second.clear();
        detail::shrink_to_fit_if_enabled(good_reads);
        detail::shrink_to_fit_if_enabled(bad_reads);
        
        good_read_map.emplace(sample_reads.first, std::move(good_reads));
        bad_read_map.emplace(std::move(sample_reads.first), std::move(bad_reads));
    }
    
    return {good_read_map, bad_read_map};
}

template <typename InputIterator>
std::vector<unsigned> positional_coverage(InputIterator first, InputIterator last, const GenomicRegion& a_region)
{
    auto num_positions = size(a_region);
    
    std::vector<unsigned> result(num_positions, 0);
    
    auto first_position = get_begin(a_region);
    
    std::for_each(first, last, [&result, first_position, num_positions] (const auto& read) {
        auto first = std::next(result.begin(), (get_begin(read) <= first_position) ? 0 : get_begin(read) - first_position);
        auto last  = std::next(result.begin(), std::min(get_end(read) - first_position, num_positions));
        std::transform(first, last, first, [] (unsigned count) { return count + 1; });
    });
    
    return result;
}

std::vector<unsigned> positional_coverage(const MappableSet<AlignedRead>& reads, const GenomicRegion& region);

template <typename T>
std::vector<unsigned> positional_coverage(const T& reads, const GenomicRegion& region)
{
    auto overlapped = overlap_range(reads.cbegin(), reads.cend(), region);
    return positional_coverage(overlapped.begin(), overlapped.end(), region);
}

template <typename T>
unsigned min_coverage(const T& reads, const GenomicRegion& region)
{
    auto positions_coverage = positional_coverage(reads, region);
    return *std::min_element(std::cbegin(positions_coverage), std::cend(positions_coverage));
}

template <typename T>
unsigned max_coverage(const T& reads, const GenomicRegion& region)
{
    auto positions_coverage = positional_coverage(reads, region);
    return *std::max_element(std::cbegin(positions_coverage), std::cend(positions_coverage));
}

template <typename T>
double mean_coverage(const T& reads, const GenomicRegion& region)
{
    auto positions_coverage = positional_coverage(reads, region);
    return mean(std::cbegin(positions_coverage), std::cend(positions_coverage));
}

template <typename T>
double stdev_coverage(const T& reads, const GenomicRegion& region)
{
    auto positions_coverage = positional_coverage(reads, region);
    return stdev(std::cbegin(positions_coverage), std::cend(positions_coverage));
}

namespace detail {
    template <typename T, typename F>
    std::unordered_map<AlignedRead, unsigned>
    f_coverage_in_read_region(const T& reads, const GenomicRegion& region, F f) {
        std::unordered_map<AlignedRead, unsigned> result {};
        result.reserve(reads.size());
        
        auto position_coverages = positional_coverage(reads, region);
        auto first_position     = get_begin(region);
        auto num_positions      = size(region);
        
        for (const auto& read : reads) {
            auto first = std::next(std::cbegin(position_coverages), (get_begin(read) <= first_position) ? 0 : get_begin(read) - first_position);
            auto last  = std::next(std::cbegin(position_coverages), std::min(get_end(read) - first_position, num_positions));
            result.emplace(read, *f(first, last));
        }
        
        return result;
    }
}

template <typename T>
std::unordered_map<AlignedRead, unsigned> min_coverage_in_read_region(const T& reads, const GenomicRegion& region)
{
    return detail::f_coverage_in_read_region(reads, region, std::min_element);
}

template <typename T>
std::unordered_map<AlignedRead, unsigned> max_coverage_in_read_region(const T& reads, const GenomicRegion& region)
{
    return detail::f_coverage_in_read_region(reads, region, std::max_element);
}

template <typename T, typename Container>
unsigned min_coverage(const ReadMap<T, Container>& reads, const GenomicRegion& region)
{
    std::vector<unsigned> sample_min_coverages(reads.size(), 0);
    
    std::transform(std::cbegin(reads), std::cend(reads), sample_min_coverages.begin(),
                   [&region] (const auto& sample_reads) {
                       return min_coverage(sample_reads.second, region);
                   });
    
    return *std::min_element(sample_min_coverages.cbegin(), sample_min_coverages.cend());
}

template <typename T, typename Container>
unsigned max_coverage(const ReadMap<T, Container>& reads, const GenomicRegion& region)
{
    std::vector<unsigned> sample_max_coverages(reads.size(), 0);
    
    std::transform(std::cbegin(reads), std::cend(reads), sample_max_coverages.begin(),
                   [&region] (const auto& sample_reads) {
                       return max_coverage(sample_reads.second, region);
                   });
    
    return *std::max_element(sample_max_coverages.cbegin(), sample_max_coverages.cend());
}


template <typename T>
std::vector<GenomicRegion> find_high_coverage_regions(const T& reads,
                                                      const GenomicRegion& region, unsigned maximum_coverage)
{
    std::vector<GenomicRegion> result {};
    
    auto positions_coverage = positional_coverage(reads, region);
    
    using Iterator = typename decltype(positions_coverage)::const_iterator;
    Iterator first {positions_coverage.cbegin()};
    Iterator current {first};
    Iterator last {positions_coverage.cend()};
    Iterator high_range_first, high_range_last;
    GenomicRegion::SizeType high_range_begin, high_range_end;
    
    auto f_is_high_coverage = [maximum_coverage] (unsigned coverage) {
        return coverage > maximum_coverage;
    };
    
    while (current != last) {
        high_range_first = std::find_if(current, last, f_is_high_coverage);
        
        if (high_range_first == last) break;
        
        high_range_last = std::find_if_not(high_range_first, last, f_is_high_coverage);
        
        high_range_begin = get_begin(region) + static_cast<GenomicRegion::SizeType>(std::distance(first, high_range_first));
        high_range_end   = high_range_begin + static_cast<GenomicRegion::SizeType>(std::distance(high_range_first, high_range_last));
        
        result.emplace_back(get_contig_name(region), high_range_begin, high_range_end);
        
        current = high_range_last;
    }
    
    result.shrink_to_fit();
    
    return result;
}

template <typename SampleIdType, typename T>
std::unordered_map<SampleIdType, std::vector<GenomicRegion>>
find_high_coverage_regions(const ReadMap<SampleIdType, T>& reads,
                           const GenomicRegion& region, unsigned maximum_coverage)
{
    std::unordered_map<SampleIdType, std::vector<GenomicRegion>> result {};
    result.reserve(reads.size());
    
    for (const auto& sample_reads : reads) {
        result.emplace(sample_reads.first, find_high_coverage_regions(sample_reads.second, region, maximum_coverage));
    }
    
    return result;
}

template <typename T>
std::vector<GenomicRegion>
find_good_coverage_regions_containing_high_coverage_positions(const T& reads, const GenomicRegion& region,
                                                              unsigned maximum_coverage, unsigned minimum_downsample_coverage)
{
    auto above_max_coverage_regions = find_high_coverage_regions(reads, region, maximum_coverage);
    
    std::vector<GenomicRegion> result {};
    
    if (above_max_coverage_regions.empty()) return result;
    
    result.reserve(above_max_coverage_regions.size());
    
    auto above_min_coverage_regions = find_high_coverage_regions(reads, region, minimum_downsample_coverage);
    
    std::copy_if(std::cbegin(above_min_coverage_regions), std::cend(above_min_coverage_regions), std::back_inserter(result),
                 [&above_max_coverage_regions] (const auto& r) {
                     return has_contained(std::cbegin(above_max_coverage_regions), std::cend(above_max_coverage_regions), r);
                 });
    
    result.shrink_to_fit();
    
    return result;
}

MappableSet<AlignedRead> downsample(const MappableSet<AlignedRead>& reads, unsigned maximum_coverage,
                                    unsigned minimum_downsample_coverage);

//template <typename T, typename Container>
//ReadMap<T, Container> downsample(ReadMap<T, Container>&& reads, unsigned max_coverage_per_sample)
//{
//    return reads;
//}

#endif /* defined(__Octopus__read_utils__) */
