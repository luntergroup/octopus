//
//  read_utils.hpp
//  Octopus
//
//  Created by Daniel Cooke on 09/04/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__read_utils__
#define __Octopus__read_utils__

#include <vector>
#include <unordered_map>
#include <iterator>  // std::begin, std::end, std::make_move_iterator, std::back_inserter, std::distance
#include <utility>   // std::move
#include <algorithm> // std::min_element, std::max_element, std::transform, std::for_each, std::min,
                     // std::find_if, std::find_if_not, std::copy_if, std::any_of, std::count_if
#include <numeric>   // std::accumulate

#include "aligned_read.hpp"
#include "read_filter.hpp"
#include "read_transform.hpp"
#include "context_iterators.hpp"
#include "mappable_algorithms.hpp"
#include "maths.hpp"
#include "mappable_set.hpp"
#include "mappable_map.hpp"

namespace Octopus
{

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
    
    filter.filter_reads(std::cbegin(reads), std::cend(reads), ContextInserter(good_reads),
                        ContextInserter(bad_reads));
    
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
    
    filter.filter_reads(std::make_move_iterator(std::begin(reads)),
                        std::make_move_iterator(std::end(reads)),
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
    
    for (auto&& sample_reads : reads) {
        auto sample_filtered_reads = filter_reads(std::move(sample_reads.second), filter);
        good_reads.emplace(sample_reads.first, std::move(sample_filtered_reads.first));
        bad_reads.emplace(sample_reads.first, std::move(sample_filtered_reads.second));
    }
    
    return {good_reads, bad_reads};
}

template <typename T, typename Container, typename ReadFilter>
auto
filter_reads(std::unordered_map<T, Container>&& reads, ReadFilter& read_filter)
{
    std::unordered_map<T, Container> good_read_map {}, bad_read_map {};
    good_read_map.reserve(reads.size());
    bad_read_map.reserve(reads.size());
    
    for (auto&& sample_reads : reads) {
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
    
    return std::make_pair(good_read_map, bad_read_map);
}

template <typename ReadMap>
void transform_reads(ReadMap& reads, ReadTransform& transformer)
{
    for (auto& p : reads) transformer.transform_reads(std::begin(p.second), std::end(p.second));
}

template <typename InputIterator>
std::vector<unsigned> positional_coverage(InputIterator first, InputIterator last, const GenomicRegion& region)
{
    auto num_positions = size(region);
    
    std::vector<unsigned> result(num_positions, 0);
    
    auto first_position = get_begin(region);
    
    std::for_each(first, last, [&result, first_position, num_positions] (const auto& read) {
        auto first = std::next(std::begin(result), (get_begin(read) <= first_position) ? 0 : get_begin(read) - first_position);
        auto last  = std::next(std::begin(result), std::min(get_end(read) - first_position, num_positions));
        std::transform(first, last, first, [] (auto count) { return count + 1; });
    });
    
    return result;
}

inline std::vector<unsigned> positional_coverage(const MappableSet<AlignedRead>& reads, const GenomicRegion& region)
{
    auto overlapped = reads.overlap_range(region);
    return positional_coverage(overlapped.begin(), overlapped.end(), region);
}

template <typename T>
std::vector<unsigned> positional_coverage(const T& reads, const GenomicRegion& region)
{
    auto overlapped = overlap_range(std::cbegin(reads), std::cend(reads), region);
    return positional_coverage(overlapped.begin(), overlapped.end(), region);
}

namespace detail
{
    template <bool> struct IsMapType {};
    using MapTypeTag    = IsMapType<true>;
    using NonMapTypeTag = IsMapType<false>;
    
    template <typename T>
    bool has_coverage(const T& reads, const GenomicRegion& region, NonMapTypeTag)
    {
        auto positions_coverage = positional_coverage(reads, region);
        return std::any_of(std::cbegin(positions_coverage), std::cend(positions_coverage),
                           [] (auto coverage) { return coverage > 0; });
    }
    
    template <typename T>
    unsigned min_coverage(const T& reads, const GenomicRegion& region, NonMapTypeTag)
    {
        if (reads.empty() || empty(region)) return 0;
        auto positions_coverage = positional_coverage(reads, region);
        return *std::min_element(std::cbegin(positions_coverage), std::cend(positions_coverage));
    }
    
    template <typename T>
    unsigned max_coverage(const T& reads, const GenomicRegion& region, NonMapTypeTag)
    {
        if (reads.empty() || empty(region)) return 0;
        auto positions_coverage = positional_coverage(reads, region);
        return *std::max_element(std::cbegin(positions_coverage), std::cend(positions_coverage));
    }
    
    template <typename T>
    double mean_coverage(const T& reads, const GenomicRegion& region, NonMapTypeTag)
    {
        if (reads.empty() || empty(region)) return 0;
        auto positions_coverage = positional_coverage(reads, region);
        return Maths::mean(positions_coverage);
    }
    
    template <typename T>
    double stdev_coverage(const T& reads, const GenomicRegion& region, NonMapTypeTag)
    {
        if (reads.empty() || empty(region)) return 0;
        auto positions_coverage = positional_coverage(reads, region);
        return Maths::stdev(positions_coverage);
    }
    
    template <typename T>
    size_t count_base_pairs(const T& reads, NonMapTypeTag)
    {
        return std::accumulate(std::cbegin(reads), std::cend(reads), size_t {},
                               [] (auto curr, const auto& read) {
                                   return curr + read.get_sequence_size();
                               });
    }
    
    template <typename T>
    size_t count_base_pairs(const T& reads, const GenomicRegion& region, NonMapTypeTag)
    {
        auto overlapped = overlap_range(std::cbegin(reads), std::cend(reads), region);
        return std::accumulate(std::cbegin(overlapped), std::cend(overlapped), size_t {},
                               [] (auto curr, const auto& read) {
                                   return curr + read.get_sequence_size();
                               });
    }
    
    template <typename T>
    size_t count_reads(const T& reads, NonMapTypeTag)
    {
        return reads.size();
    }
    
    template <typename T>
    size_t count_reads(const T& reads, const GenomicRegion& region, NonMapTypeTag)
    {
        return count_overlapped(std::cbegin(reads), std::cend(reads), region);
    }
    
    inline size_t count_reads(const MappableSet<AlignedRead>& reads, const GenomicRegion& region, NonMapTypeTag)
    {
        return reads.count_overlapped(region);
    }
    
    template <typename T>
    size_t count_forward(const T& reads, NonMapTypeTag)
    {
        return std::count_if(std::cbegin(reads), std::cend(reads),
                             [] (const auto& read) { return !read.is_marked_reverse_mapped(); });
    }
    
    
    template <typename T>
    size_t count_forward(const T& reads, const GenomicRegion& region, NonMapTypeTag)
    {
        auto overlapped = overlap_range(std::cbegin(reads), std::cend(reads), region);
        return std::count_if(std::cbegin(overlapped), std::cend(overlapped),
                             [] (const auto& read) { return !read.is_marked_reverse_mapped(); });
    }
    
    inline size_t count_forward(const MappableSet<AlignedRead>& reads, const GenomicRegion& region, NonMapTypeTag)
    {
        auto overlapped = reads.overlap_range(region);
        return std::count_if(std::cbegin(overlapped), std::cend(overlapped),
                             [] (const auto& read) { return !read.is_marked_reverse_mapped(); });
    }
    
    template <typename T>
    size_t count_reverse(const T& reads, NonMapTypeTag)
    {
        return std::count_if(std::cbegin(reads), std::cend(reads),
                             [] (const auto& read) { return read.is_marked_reverse_mapped(); });
    }
    
    template <typename T>
    size_t count_reverse(const T& reads, const GenomicRegion& region, NonMapTypeTag)
    {
        auto overlapped = overlap_range(std::cbegin(reads), std::cend(reads), region);
        return std::count_if(std::cbegin(overlapped), std::cend(overlapped),
                               [] (const auto& read) { return read.is_marked_reverse_mapped(); });
    }
    
    inline size_t count_reverse(const MappableSet<AlignedRead>& reads, const GenomicRegion& region, NonMapTypeTag)
    {
        auto overlapped = reads.overlap_range(region);
        return std::count_if(std::cbegin(overlapped), std::cend(overlapped),
                             [] (const auto& read) { return read.is_marked_reverse_mapped(); });
    }
    
    template <typename T>
    size_t count_mapq_zero(const T& reads, NonMapTypeTag)
    {
        return std::count_if(std::cbegin(reads), std::cend(reads),
                             [] (const auto& read) { return read.get_mapping_quality() == 0; });
    }
    
    template <typename T>
    size_t count_mapq_zero(const T& reads, const GenomicRegion& region, NonMapTypeTag)
    {
        auto overlapped = overlap_range(std::cbegin(reads), std::cend(reads), region);
        return std::count_if(std::cbegin(overlapped), std::cend(overlapped),
                             [] (const auto& read) { return read.get_mapping_quality() == 0; });
    }
    
    inline size_t count_mapq_zero(const MappableSet<AlignedRead>& reads, const GenomicRegion& region, NonMapTypeTag)
    {
        auto overlapped = reads.overlap_range(region);
        return std::count_if(std::cbegin(overlapped), std::cend(overlapped),
                             [] (const auto& read) { return read.get_mapping_quality() == 0; });
    }
    
    template <typename T>
    double rmq_mapping_quality(const T& reads, NonMapTypeTag)
    {
        std::vector<double> qualities(reads.size());
        
        std::transform(std::cbegin(reads), std::cend(reads), std::begin(qualities),
                       [] (const auto& read) { return static_cast<double>(read.get_mapping_quality()); });
        
        return Maths::rmq<double>(qualities);
    }
    
    template <typename T>
    double rmq_mapping_quality(const T& reads, const GenomicRegion& region, NonMapTypeTag)
    {
        auto overlapped = overlap_range(std::cbegin(reads), std::cend(reads), region);
        std::vector<double> qualities(size(overlapped));
        
        std::transform(std::cbegin(overlapped), std::cend(overlapped), std::begin(qualities),
                       [] (const auto& read) { return static_cast<double>(read.get_mapping_quality()); });
        
        return Maths::rmq<double>(qualities);
    }
    
    inline double rmq_mapping_quality(const MappableSet<AlignedRead>& reads, const GenomicRegion& region, NonMapTypeTag)
    {
        auto overlapped = reads.overlap_range(region);
        std::vector<double> qualities(size(overlapped));
        
        std::transform(std::cbegin(overlapped), std::cend(overlapped), std::begin(qualities),
                       [] (const auto& read) { return static_cast<double>(read.get_mapping_quality()); });
        
        return Maths::rmq<double>(qualities);
    }
    
    template <typename T>
    double rmq_base_quality(const T& reads, NonMapTypeTag)
    {
        std::vector<double> qualities {};
        qualities.reserve(count_base_pairs(reads, NonMapTypeTag()));
        
        for (const auto& read : reads) {
            for (const auto quality : read.get_qualities()) {
                qualities.push_back(static_cast<double>(quality));
            }
        }
        
        return Maths::rmq<double>(qualities);
    }
    
    template <typename T>
    double rmq_base_quality(const T& reads, const GenomicRegion& region, NonMapTypeTag)
    {
        auto overlapped = overlap_range(std::cbegin(reads), std::cend(reads), region);
        
        std::vector<double> qualities {};
        qualities.reserve(count_base_pairs(reads, region, NonMapTypeTag()));
        
        std::for_each(std::cbegin(overlapped), std::cend(overlapped), [&qualities] (const auto& read) {
            for (const auto quality : read.get_qualities()) {
                qualities.push_back(static_cast<double>(quality));
            }
        });
        
        return Maths::rmq<double>(qualities);
    }
    
    inline double rmq_base_quality(const MappableSet<AlignedRead>& reads, const GenomicRegion& region, NonMapTypeTag)
    {
        auto overlapped = reads.overlap_range(region);
        
        std::vector<double> qualities {};
        qualities.reserve(count_base_pairs(reads, region, NonMapTypeTag()));
        
        std::for_each(std::cbegin(overlapped), std::cend(overlapped), [&qualities] (const auto& read) {
            for (const auto quality : read.get_qualities()) {
                qualities.push_back(static_cast<double>(quality));
            }
        });
        
        return Maths::rmq<double>(qualities);
    }
    
    template <typename T>
    size_t count_base_pairs(const T& reads, MapTypeTag)
    {
        return std::accumulate(std::cbegin(reads), std::cend(reads), size_t {},
                               [] (auto curr, const auto& sample_reads) {
                                   return curr + count_base_pairs(sample_reads.second, NonMapTypeTag());
                               });
    }
    
    template <typename T>
    size_t count_base_pairs(const T& reads, const GenomicRegion& region, MapTypeTag)
    {
        return std::accumulate(std::cbegin(reads), std::cend(reads), size_t {},
                               [&region] (auto curr, const auto& sample_reads) {
                                   return curr + count_base_pairs(sample_reads.second, region, NonMapTypeTag());
                               });
    }
    
    template <typename ReadMap>
    bool has_coverage(const ReadMap& reads, MapTypeTag)
    {
        return std::any_of(std::cbegin(reads), std::cend(reads),
                           [] (const auto& sample_reads) {
                               return has_coverage(sample_reads.second, NonMapTypeTag());
                           });
    }
    
    template <typename ReadMap>
    bool has_coverage(const ReadMap& reads, const GenomicRegion& region, MapTypeTag)
    {
        return std::any_of(std::cbegin(reads), std::cend(reads),
                           [&region] (const auto& sample_reads) {
                               return has_coverage(sample_reads.second, region, NonMapTypeTag());
                           });
    }
    
    template <typename ReadMap>
    unsigned min_coverage(const ReadMap& reads, MapTypeTag)
    {
        if (reads.empty()) return 0;
        
        std::vector<unsigned> sample_min_coverages(reads.size(), 0);
        
        std::transform(std::cbegin(reads), std::cend(reads), std::begin(sample_min_coverages),
                       [] (const auto& sample_reads) {
                           return min_coverage(sample_reads.second, NonMapTypeTag());
                       });
        
        return *std::min_element(cbegin(sample_min_coverages), std::cend(sample_min_coverages));
    }
    
    template <typename ReadMap>
    unsigned min_coverage(const ReadMap& reads, const GenomicRegion& region, MapTypeTag)
    {
        if (reads.empty()) return 0;
        
        std::vector<unsigned> sample_min_coverages(reads.size(), 0);
        
        std::transform(std::cbegin(reads), std::cend(reads), std::begin(sample_min_coverages),
                       [&region] (const auto& sample_reads) {
                           return min_coverage(sample_reads.second, region, NonMapTypeTag());
                       });
        
        return *std::min_element(cbegin(sample_min_coverages), std::cend(sample_min_coverages));
    }
    
    template <typename ReadMap>
    unsigned max_coverage(const ReadMap& reads, const GenomicRegion& region, MapTypeTag)
    {
        if (reads.empty()) return 0;
        
        std::vector<unsigned> sample_max_coverages(reads.size(), 0);
        
        std::transform(std::cbegin(reads), std::cend(reads), std::begin(sample_max_coverages),
                       [&region] (const auto& sample_reads) {
                           return max_coverage(sample_reads.second, region, NonMapTypeTag());
                       });
        
        return *std::max_element(std::cbegin(sample_max_coverages), std::cend(sample_max_coverages));
    }
    
    template <typename T>
    double mean_coverage(const T& reads, const GenomicRegion& region, MapTypeTag)
    {
        if (reads.empty()) return 0;
        
        std::vector<double> sample_mean_coverages(reads.size(), 0);
        
        std::transform(std::cbegin(reads), std::cend(reads), std::begin(sample_mean_coverages),
                       [&region] (const auto& sample_reads) {
                           return mean(sample_reads.second, region, NonMapTypeTag());
                       });
        
        return Maths::mean(sample_mean_coverages);
    }
    
    template <typename T>
    double stdev_coverage(const T& reads, const GenomicRegion& region, MapTypeTag)
    {
        if (reads.empty()) return 0;
        
        std::vector<double> sample_stdev_coverages(reads.size(), 0);
        
        std::transform(std::cbegin(reads), std::cend(reads), std::begin(sample_stdev_coverages),
                       [&region] (const auto& sample_reads) {
                           return stdev(sample_reads.second, region, NonMapTypeTag());
                       });
        
        return Maths::stdev(sample_stdev_coverages);
    }
    
    template <typename T>
    size_t count_reads(const T& reads, MapTypeTag)
    {
        return Maths::sum_sizes(reads);
    }
    
    template <typename T>
    size_t count_reads(const T& reads, const GenomicRegion& region, MapTypeTag)
    {
        return std::accumulate(std::cbegin(reads), std::cend(reads), 0,
                               [&region] (auto curr, const auto& sample_reads) {
                                   return curr + count_reads(sample_reads.second, region, NonMapTypeTag());
                               });
    }
    
    template <typename T>
    size_t count_forward(const T& reads, MapTypeTag)
    {
        return std::accumulate(std::cbegin(reads), std::cend(reads), 0,
                               [] (auto curr, const auto& sample_reads) {
                                   return curr + count_forward(sample_reads.second, NonMapTypeTag());
                               });
    }
    
    template <typename T>
    size_t count_forward(const T& reads, const GenomicRegion& region, MapTypeTag)
    {
        return std::accumulate(std::cbegin(reads), std::cend(reads), 0,
                               [&region] (auto curr, const auto& sample_reads) {
                                   return curr + count_forward(sample_reads.second, region, NonMapTypeTag());
                               });
    }
    
    template <typename T>
    size_t count_reverse(const T& reads, MapTypeTag)
    {
        return std::accumulate(std::cbegin(reads), std::cend(reads), 0,
                               [] (auto curr, const auto& sample_reads) {
                                   return curr + count_reverse(sample_reads.second, NonMapTypeTag());
                               });
    }
    
    template <typename T>
    size_t count_reverse(const T& reads, const GenomicRegion& region, MapTypeTag)
    {
        return std::accumulate(std::cbegin(reads), std::cend(reads), 0,
                               [&region] (auto curr, const auto& sample_reads) {
                                   return curr + count_reverse(sample_reads.second, region, NonMapTypeTag());
                               });
    }
    
    template <typename T>
    size_t count_mapq_zero(const T& reads, MapTypeTag)
    {
        return std::accumulate(std::cbegin(reads), std::cend(reads), 0,
                               [] (auto curr, const auto& sample_reads) {
                                   return curr + count_mapq_zero(sample_reads.second, NonMapTypeTag());
                               });
    }
    
    template <typename T>
    size_t count_mapq_zero(const T& reads, const GenomicRegion& region, MapTypeTag)
    {
        return std::accumulate(std::cbegin(reads), std::cend(reads), 0,
                               [&region] (auto curr, const auto& sample_reads) {
                                   return curr + count_mapq_zero(sample_reads.second, region, NonMapTypeTag());
                               });
    }
    
    template <typename T>
    double rmq_mapping_quality(const T& reads, MapTypeTag)
    {
        std::vector<double> qualities {};
        qualities.reserve(Maths::sum_sizes(reads));
        
        for (const auto& sample_reads : reads) {
            std::transform(std::cbegin(sample_reads.second), std::cend(sample_reads.second),
                           std::back_inserter(qualities), [] (const auto& read) {
                               return static_cast<double>(read.get_mapping_quality());
                           });
        }
        
        return Maths::rmq<double>(qualities);
    }
    
    template <typename T>
    double rmq_mapping_quality(const T& reads, const GenomicRegion& region, MapTypeTag)
    {
        std::vector<double> qualities {};
        qualities.reserve(Maths::sum_sizes(reads));
        
        for (const auto& sample_reads : reads) {
            auto overlapped = overlap_range(std::cbegin(sample_reads.second),
                                            std::cend(sample_reads.second), region);
            std::transform(std::cbegin(overlapped), std::cend(overlapped),
                           std::back_inserter(qualities), [] (const auto& read) {
                               return static_cast<double>(read.get_mapping_quality());
                           });
        }
        
        return Maths::rmq<double>(qualities);
    }
    
    template <typename T>
    double rmq_base_quality(const T& reads, MapTypeTag)
    {
        std::vector<double> qualities {};
        qualities.reserve(count_base_pairs(reads, MapTypeTag()));
        
        for (const auto& sample_reads : reads) {
            for (const auto& read : sample_reads.second) {
                for (const auto quality : read.get_qualities()) {
                    qualities.push_back(static_cast<double>(quality));
                }
            }
        }
        
        return Maths::rmq<double>(qualities);
    }
    
    template <typename T>
    double rmq_base_quality(const T& reads, const GenomicRegion& region, MapTypeTag)
    {
        std::vector<double> qualities {};
        qualities.reserve(count_base_pairs(reads, region, MapTypeTag()));
        
        for (const auto& sample_reads : reads) {
            auto overlapped = overlap_range(std::cbegin(sample_reads.second),
                                            std::cend(sample_reads.second), region);
            std::for_each(std::cbegin(overlapped), std::cend(overlapped), [&qualities] (const auto& read) {
                for (const auto quality : read.get_qualities()) {
                    qualities.push_back(static_cast<double>(quality));
                }
            });
        }
        
        return Maths::rmq<double>(qualities);
    }
} // namespace detail

template <typename T>
bool has_coverage(const T& reads, const GenomicRegion& region)
{
    return detail::has_coverage(reads, region, detail::IsMapType<!std::is_same<typename T::value_type, AlignedRead>::value>());
}

template <typename T>
unsigned min_coverage(const T& reads, const GenomicRegion& region)
{
    return detail::min_coverage(reads, region, detail::IsMapType<!std::is_same<typename T::value_type, AlignedRead>::value>());
}

template <typename T>
unsigned max_coverage(const T& reads, const GenomicRegion& region)
{
    return detail::max_coverage(reads, region, detail::IsMapType<!std::is_same<typename T::value_type, AlignedRead>::value>());
}

template <typename T>
double mean_coverage(const T& reads, const GenomicRegion& region)
{
    return detail::mean_coverage(reads, region, detail::IsMapType<!std::is_same<typename T::value_type, AlignedRead>::value>());
}

template <typename T>
double stdev_coverage(const T& reads, const GenomicRegion& region)
{
    return detail::stdev_coverage(reads, region, detail::IsMapType<!std::is_same<typename T::value_type, AlignedRead>::value>());
}

template <typename ReadMap>
size_t count_samples_with_coverage(const ReadMap& reads, const GenomicRegion& region)
{
    return std::count_if(std::cbegin(reads), std::cend(reads),
                         [&region] (const auto& sample_reads) {
                             return has_coverage(sample_reads.second, region);
                         });
}

template <typename ReadMap>
unsigned sum_max_coverages(const ReadMap& reads, const GenomicRegion& region)
{
    return std::accumulate(std::cbegin(reads), std::cend(reads), 0,
                           [&region] (auto curr, const auto& sample_reads) {
                               return curr + max_coverage(sample_reads.second, region);
                           });
}

template <typename T>
bool count_base_pairs(const T& reads)
{
    return detail::count_base_pairs(reads, detail::IsMapType<!std::is_same<typename T::value_type, AlignedRead>::value>());
}

template <typename T>
bool count_base_pairs(const T& reads, const GenomicRegion& region)
{
    return detail::count_base_pairs(reads, region, detail::IsMapType<!std::is_same<typename T::value_type, AlignedRead>::value>());
}

template <typename T>
size_t count_reads(const T& reads)
{
    return detail::count_reads(reads, detail::IsMapType<!std::is_same<typename T::value_type, AlignedRead>::value>());
}

template <typename T>
size_t count_reads(const T& reads, const GenomicRegion& region)
{
    return detail::count_reads(reads, region, detail::IsMapType<!std::is_same<typename T::value_type, AlignedRead>::value>());
}

template <typename T>
size_t count_forward(const T& reads)
{
    return detail::count_forward(reads, detail::IsMapType<!std::is_same<typename T::value_type, AlignedRead>::value>());
}

template <typename T>
size_t count_forward(const T& reads, const GenomicRegion& region)
{
    return detail::count_forward(reads, region, detail::IsMapType<!std::is_same<typename T::value_type, AlignedRead>::value>());
}

template <typename T>
size_t count_reverse(const T& reads)
{
    return detail::count_reverse(reads, detail::IsMapType<!std::is_same<typename T::value_type, AlignedRead>::value>());
}

template <typename T>
size_t count_reverse(const T& reads, const GenomicRegion& region)
{
    return detail::count_reverse(reads, region, detail::IsMapType<!std::is_same<typename T::value_type, AlignedRead>::value>());
}

template <typename T>
double strand_bias(const T& reads)
{
    auto num_forward_reads = count_forward(reads);
    auto num_reverse_reads = count_reverse(reads);
    return static_cast<double>(num_forward_reads) / num_reverse_reads;
}

template <typename T>
double strand_bias(const T& reads, const GenomicRegion& region)
{
    auto num_forward_reads = count_forward(reads, region);
    auto num_reverse_reads = count_reverse(reads, region);
    return static_cast<double>(num_forward_reads) / num_reverse_reads;
}

template <typename T>
size_t count_mapq_zero(const T& reads)
{
    return detail::count_mapq_zero(reads, detail::IsMapType<!std::is_same<typename T::value_type, AlignedRead>::value>());
}

template <typename T>
size_t count_mapq_zero(const T& reads, const GenomicRegion& region)
{
    return detail::count_mapq_zero(reads, region, detail::IsMapType<!std::is_same<typename T::value_type, AlignedRead>::value>());
}

template <typename T>
double rmq_mapping_quality(const T& reads)
{
    return detail::rmq_mapping_quality(reads, detail::IsMapType<!std::is_same<typename T::value_type, AlignedRead>::value>());
}

template <typename T>
double rmq_mapping_quality(const T& reads, const GenomicRegion& region)
{
    return detail::rmq_mapping_quality(reads, region, detail::IsMapType<!std::is_same<typename T::value_type, AlignedRead>::value>());
}

template <typename T>
double rmq_base_quality(const T& reads)
{
    return detail::rmq_base_quality(reads, detail::IsMapType<!std::is_same<typename T::value_type, AlignedRead>::value>());
}

template <typename T>
double rmq_base_quality(const T& reads, const GenomicRegion& region)
{
    return detail::rmq_base_quality(reads, region, detail::IsMapType<!std::is_same<typename T::value_type, AlignedRead>::value>());
}

namespace detail
{
    template <typename T, typename F>
    std::unordered_map<AlignedRead, unsigned>
    coverages_in_read_regions(const T& reads, const GenomicRegion& region, F f) {
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
} // namespace detail

template <typename T>
std::unordered_map<AlignedRead, unsigned> get_min_coverages_in_read_regions(const T& reads, const GenomicRegion& region)
{
    return detail::coverages_in_read_regions(reads, region, std::min_element);
}

template <typename T>
std::unordered_map<AlignedRead, unsigned> get_max_coverages_in_read_regions(const T& reads, const GenomicRegion& region)
{
    return detail::coverages_in_read_regions(reads, region, std::max_element);
}

template <typename T>
std::vector<GenomicRegion>
find_high_coverage_regions(const T& reads, const GenomicRegion& region, const unsigned max_coverage)
{
    using SizeType = GenomicRegion::SizeType;
    
    std::vector<GenomicRegion> result {};
    
    auto positions_coverage = positional_coverage(reads, region);
    
    using Iterator = typename decltype(positions_coverage)::const_iterator;
    
    Iterator first {positions_coverage.cbegin()};
    Iterator current {first};
    Iterator last {positions_coverage.cend()};
    Iterator high_range_first, high_range_last;
    
    SizeType high_range_begin, high_range_end;
    
    while (current != last) {
        auto is_high_coverage = [max_coverage] (unsigned coverage) { return coverage > max_coverage; };
        
        high_range_first = std::find_if(current, last, is_high_coverage);
        
        if (high_range_first == last) break;
        
        high_range_last = std::find_if_not(high_range_first, last, is_high_coverage);
        
        high_range_begin = get_begin(region) + static_cast<SizeType>(std::distance(first, high_range_first));
        high_range_end   = high_range_begin + static_cast<SizeType>(std::distance(high_range_first, high_range_last));
        
        result.emplace_back(get_contig_name(region), high_range_begin, high_range_end);
        
        current = high_range_last;
    }
    
    result.shrink_to_fit();
    
    return result;
}

template <typename ReadMap>
std::unordered_map<typename ReadMap::key_type, std::vector<GenomicRegion>>
find_high_coverage_regions(const ReadMap& reads, const GenomicRegion& region, const unsigned max_coverage)
{
    std::unordered_map<typename ReadMap::key_type, std::vector<GenomicRegion>> result {};
    result.reserve(reads.size());
    
    for (const auto& sample_reads : reads) {
        result.emplace(sample_reads.first, find_high_coverage_regions(sample_reads.second, region, max_coverage));
    }
    
    return result;
}

MappableSet<AlignedRead>
downsample(const MappableSet<AlignedRead>& reads, const unsigned max_coverage, const unsigned min_downsample_coverage);

template <typename T>
MappableMap<T, AlignedRead>
downsample(const MappableMap<T, AlignedRead>& reads, const unsigned max_coverage, const unsigned min_downsample_coverage)
{
    MappableMap<T, AlignedRead> result {};
    result.reserve(reads.size());
    
    for (const auto& sample_reads : reads) {
        result.emplace(sample_reads.first, downsample(sample_reads.second, max_coverage, min_downsample_coverage));
    }
    
    return result;
}

// TODO
AlignedRead find_next_segment(const AlignedRead& read, const MappableMap<GenomicRegion::StringType, AlignedRead>& reads);
    
// TODO
MappableSet<AlignedRead> find_chimeras(const AlignedRead& read, const MappableSet<AlignedRead>& reads);

} // namespace Octopus

#endif /* defined(__Octopus__read_utils__) */
