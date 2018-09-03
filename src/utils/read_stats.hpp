// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef read_stats_hpp
#define read_stats_hpp

#include <vector>
#include <iterator>
#include <utility>
#include <algorithm>
#include <numeric>
#include <type_traits>

#include <boost/range/algorithm.hpp>
#include <boost/optional.hpp>

#include "basics/aligned_read.hpp"
#include "containers/mappable_map.hpp"
#include "mappable_algorithms.hpp"
#include "maths.hpp"
#include "type_tricks.hpp"

namespace octopus {

namespace detail {

template <typename T>
constexpr bool is_aligned_read = std::is_same<std::decay_t<T>, AlignedRead>::value;

template <typename Container>
constexpr bool is_aligned_read_container = is_aligned_read<typename Container::value_type>;

struct IsForward
{
    bool operator()(const AlignedRead& read) const
    {
        return is_forward_strand(read);
    }
};

struct IsReverse
{
    bool operator()(const AlignedRead& read) const
    {
        return is_reverse_strand(read);
    }
};

struct IsMappingQualityZero
{
    bool operator()(const AlignedRead& read) const noexcept
    {
        return read.mapping_quality() == 0;
    }
};

template <typename Map>
bool all_empty(const Map& map) noexcept
{
    return std::all_of(std::cbegin(map), std::cend(map), [] (const auto& p) { return p.second.empty(); });
}

template <typename T>
bool has_coverage(const T& reads, NonMapTag)
{
    using octopus::has_coverage;
    return has_coverage(reads);
}

template <typename T>
bool has_coverage(const T& reads, const GenomicRegion& region, NonMapTag)
{
    using octopus::has_coverage;
    return has_coverage(reads, region);
}

template <typename T>
unsigned min_coverage(const T& reads, NonMapTag)
{
    using octopus::min_coverage;
    return min_coverage(reads);
}

template <typename T>
unsigned min_coverage(const T& reads, const GenomicRegion& region, NonMapTag)
{
    using octopus::min_coverage;
    return min_coverage(reads, region);
}

template <typename T>
unsigned max_coverage(const T& reads, NonMapTag)
{
    using octopus::max_coverage;
    return max_coverage(reads);
}

template <typename T>
unsigned max_coverage(const T& reads, const GenomicRegion& region, NonMapTag)
{
    using octopus::max_coverage;
    return max_coverage(reads, region);
}

template <typename T>
double mean_coverage(const T& reads, const GenomicRegion& region, NonMapTag)
{
    if (reads.empty() || is_empty(region)) return 0;
    return maths::mean(calculate_positional_coverage(reads, region));
}

template <typename T>
double mean_coverage(const T& reads, NonMapTag)
{
    if (reads.empty()) return 0;
    return mean_coverage(reads, encompassing_region(reads), NonMapTag {});
}

template <typename T>
double stdev_coverage(const T& reads, const GenomicRegion& region, NonMapTag)
{
    if (reads.empty() || is_empty(region)) return 0;
    return maths::stdev(calculate_positional_coverage(reads, region));
}

template <typename T>
double stdev_coverage(const T& reads, NonMapTag)
{
    if (reads.empty()) return 0;
    return stdev_coverage(reads, encompassing_region(reads), NonMapTag {});
}

template <typename T>
std::size_t count_reads(const T& reads, NonMapTag)
{
    return reads.size();
}

template <typename T>
std::size_t count_reads(const T& reads, const GenomicRegion& region, NonMapTag)
{
    return count_overlapped(reads, region);
}

template <typename T>
std::size_t count_forward(const T& reads, NonMapTag)
{
    static_assert(is_aligned_read_container<T>, "T must be a container of AlignedReads");
    return std::count_if(std::cbegin(reads), std::cend(reads), IsForward {});
}

template <typename T>
std::size_t count_forward(const T& reads, const GenomicRegion& region, NonMapTag)
{
    static_assert(is_aligned_read_container<T>, "T must be a container of AlignedReads");
    const auto overlapped = overlap_range(reads, region);
    return std::count_if(std::cbegin(overlapped), std::cend(overlapped), IsForward {});
}

template <typename T>
std::size_t count_reverse(const T& reads, NonMapTag)
{
    static_assert(is_aligned_read_container<T>, "T must be a container of AlignedReads");
    return std::count_if(std::cbegin(reads), std::cend(reads), IsReverse {});
}

template <typename T>
std::size_t count_reverse(const T& reads, const GenomicRegion& region, NonMapTag)
{
    static_assert(is_aligned_read_container<T>, "T must be a container of AlignedReads");
    const auto overlapped = overlap_range(reads, region);
    return std::count_if(std::cbegin(overlapped), std::cend(overlapped), IsReverse {});
}

template <typename T>
std::pair<std::size_t, std::size_t> count_directions(const T& reads, NonMapTag)
{
    static_assert(is_aligned_read_container<T>, "T must be a container of AlignedReads");
    const auto num_forward = std::count_if(std::cbegin(reads), std::cend(reads), IsForward {});
    return std::make_pair(num_forward, reads.size() - num_forward);
}

template <typename T>
std::pair<std::size_t, std::size_t> count_directions(const T& reads, const GenomicRegion& region, NonMapTag)
{
    static_assert(is_aligned_read_container<T>, "T must be a container of AlignedReads");
    const auto overlapped = overlap_range(reads, region);
    const auto num_forward = std::count_if(std::cbegin(overlapped), std::cend(overlapped), IsForward {});
    return std::make_pair(num_forward, size(overlapped) - num_forward);
}

struct IsShorter
{
    bool operator()(const AlignedRead& lhs, const AlignedRead& rhs) const noexcept
    {
        return sequence_size(lhs) < sequence_size(rhs);
    }
};

template <typename T>
AlignedRead::NucleotideSequence::size_type min_read_length(const T& reads, NonMapTag)
{
    if (reads.empty()) return 0;
    return sequence_size(*std::min_element(std::cbegin(reads), std::cend(reads), IsShorter {}));
}

template <typename T>
AlignedRead::NucleotideSequence::size_type min_read_length(const T& reads, const GenomicRegion& region, NonMapTag)
{
    const auto overlapped = overlap_range(reads, region);
    if (empty(overlapped)) return 0;
    return sequence_size(*std::min_element(std::cbegin(overlapped), std::cend(overlapped), IsShorter {}));
}

template <typename T>
AlignedRead::NucleotideSequence::size_type max_read_length(const T& reads, NonMapTag)
{
    if (reads.empty()) return 0;
    return sequence_size(*std::min_element(std::cbegin(reads), std::cend(reads), IsShorter {}));
}

template <typename T>
AlignedRead::NucleotideSequence::size_type max_read_length(const T& reads, const GenomicRegion& region, NonMapTag)
{
    const auto overlapped = overlap_range(reads, region);
    if (empty(overlapped)) return 0;
    return sequence_size(*std::max_element(std::cbegin(overlapped), std::cend(overlapped), IsShorter {}));
}

template <typename T>
std::size_t count_base_pairs(const T& reads, NonMapTag)
{
    static_assert(is_aligned_read_container<T>, "T must be a container of AlignedReads");
    return std::accumulate(std::cbegin(reads), std::cend(reads), std::size_t {0},
                           [] (const auto curr, const auto& read) {
                               return curr + static_cast<std::size_t>(sequence_size(read));
                           });
}

template <typename T>
std::size_t count_base_pairs(const T& reads, const GenomicRegion& region, NonMapTag)
{
    static_assert(is_aligned_read_container<T>, "T must be a container of AlignedReads");
    const auto overlapped = overlap_range(reads, region);
    return std::accumulate(std::cbegin(overlapped), std::cend(overlapped), std::size_t {0},
                           [&region] (const auto curr, const auto& read) {
                               return curr + sequence_size(read, region);
                           });
}

template <typename T>
std::size_t count_forward_base_pairs(const T& reads, NonMapTag)
{
    static_assert(is_aligned_read_container<T>, "T must be a container of AlignedReads");
    return std::accumulate(std::cbegin(reads), std::cend(reads), std::size_t {0},
                           [] (const auto curr, const auto& read) {
                               return curr + ((IsForward()(read)) ? static_cast<std::size_t>(sequence_size(read)) : 0);
                           });
}

template <typename T>
std::size_t count_forward_base_pairs(const T& reads, const GenomicRegion& region, NonMapTag)
{
    static_assert(is_aligned_read_container<T>, "T must be a container of AlignedReads");
    const auto overlapped = overlap_range(reads, region);
    return std::accumulate(std::cbegin(overlapped), std::cend(overlapped), std::size_t {0},
                           [&region] (const auto curr, const auto& read) {
                               return curr + ((IsForward()(read)) ? static_cast<std::size_t>(sequence_size(read, region)) : 0);
                           });
}

template <typename T>
std::size_t count_reverse_base_pairs(const T& reads, NonMapTag)
{
    static_assert(is_aligned_read_container<T>, "T must be a container of AlignedReads");
    return std::accumulate(std::cbegin(reads), std::cend(reads), std::size_t {},
                           [] (const auto curr, const auto& read) {
                               return curr + ((IsReverse()(read)) ? static_cast<std::size_t>(sequence_size(read)) : 0);
                           });
}

template <typename T>
std::size_t count_reverse_base_pairs(const T& reads, const GenomicRegion& region, NonMapTag)
{
    static_assert(is_aligned_read_container<T>, "T must be a container of AlignedReads");
    const auto overlapped = overlap_range(reads, region);
    return std::accumulate(std::cbegin(overlapped), std::cend(overlapped), std::size_t {0},
                           [&region] (const auto curr, const auto& read) {
                               return curr + ((IsReverse()(read)) ? static_cast<std::size_t>(sequence_size(read, region)) : 0);
                           });
}

template <typename T>
std::size_t count_mapq_zero(const T& reads, NonMapTag)
{
    static_assert(is_aligned_read_container<T>, "T must be a container of AlignedReads");
    return std::count_if(std::cbegin(reads), std::cend(reads), IsMappingQualityZero {});
}

template <typename T>
std::size_t count_mapq_zero(const T& reads, const GenomicRegion& region, NonMapTag)
{
    static_assert(is_aligned_read_container<T>, "T must be a container of AlignedReads");
    const auto overlapped = overlap_range(reads, region);
    return std::count_if(std::cbegin(overlapped), std::cend(overlapped), IsMappingQualityZero {});
}

struct MappingQualityLess
{
    bool operator()(const AlignedRead& lhs, const AlignedRead& rhs) const noexcept
    {
        return lhs.mapping_quality() < rhs.mapping_quality();
    }
};

template <typename T>
auto min_mapping_quality(const T& reads, NonMapTag)
{
    static_assert(is_aligned_read_container<T>, "T must be a container of AlignedReads");
    auto result_itr = std::min_element(std::cbegin(reads), std::cend(reads), MappingQualityLess {});
    return result_itr->mapping_quality();
}

template <typename T>
auto min_mapping_quality(const T& reads, const GenomicRegion& region, NonMapTag)
{
    static_assert(is_aligned_read_container<T>, "T must be a container of AlignedReads");
    const auto overlapped = overlap_range(reads, region);
    auto result_itr = std::min_element(std::cbegin(overlapped), std::cend(overlapped), MappingQualityLess {});
    return result_itr->mapping_quality();
}

template <typename T>
auto max_mapping_quality(const T& reads, NonMapTag)
{
    static_assert(is_aligned_read_container<T>, "T must be a container of AlignedReads");
    auto result_itr = std::max_element(std::cbegin(reads), std::cend(reads), MappingQualityLess {});
    return result_itr->mapping_quality();
}

template <typename T>
auto max_mapping_quality(const T& reads, const GenomicRegion& region, NonMapTag)
{
    static_assert(is_aligned_read_container<T>, "T must be a container of AlignedReads");
    const auto overlapped = overlap_range(reads, region);
    auto result_itr = std::max_element(std::cbegin(overlapped), std::cend(overlapped), MappingQualityLess {});
    return result_itr->mapping_quality();
}

template <typename T>
double median_mapping_quality(const T& reads, NonMapTag)
{
    static_assert(is_aligned_read_container<T>, "T must be a container of AlignedReads");
    std::vector<double> qualities(reads.size());
    std::transform(std::cbegin(reads), std::cend(reads), std::begin(qualities),
                   [] (const auto& read) { return static_cast<double>(read.mapping_quality()); });
    return maths::median<double>(qualities);
}

template <typename T>
double median_mapping_quality(const T& reads, const GenomicRegion& region, NonMapTag)
{
    static_assert(is_aligned_read_container<T>, "T must be a container of AlignedReads");
    const auto overlapped = overlap_range(reads, region);
    std::vector<double> qualities(size(overlapped));
    std::transform(std::cbegin(overlapped), std::cend(overlapped), std::begin(qualities),
                   [] (const auto& read) { return static_cast<double>(read.mapping_quality()); });
    return maths::median<double>(qualities);
}

template <typename T>
double rmq_mapping_quality(const T& reads, NonMapTag)
{
    static_assert(is_aligned_read_container<T>, "T must be a container of AlignedReads");
    std::vector<double> qualities(reads.size());
    std::transform(std::cbegin(reads), std::cend(reads), std::begin(qualities),
                   [] (const auto& read) { return static_cast<double>(read.mapping_quality()); });
    return maths::rmq<double>(qualities);
}

template <typename T>
double rmq_mapping_quality(const T& reads, const GenomicRegion& region, NonMapTag)
{
    static_assert(is_aligned_read_container<T>, "T must be a container of AlignedReads");
    const auto overlapped = overlap_range(reads, region);
    std::vector<double> qualities(size(overlapped));
    std::transform(std::cbegin(overlapped), std::cend(overlapped), std::begin(qualities),
                   [] (const auto& read) { return static_cast<double>(read.mapping_quality()); });
    return maths::rmq<double>(qualities);
}

template <typename T>
double median_base_quality(const T& reads, NonMapTag)
{
    static_assert(is_aligned_read_container<T>, "T must be a container of AlignedReads");
    std::vector<double> qualities {};
    qualities.reserve(count_base_pairs(reads, NonMapTag {}));
    for (const auto& read : reads) {
        for (const auto quality : read.base_qualities()) {
            qualities.push_back(static_cast<double>(quality));
        }
    }
    return maths::median<double>(qualities);
}

template <typename T>
double median_base_quality(const T& reads, const GenomicRegion& region, NonMapTag)
{
    static_assert(is_aligned_read_container<T>, "T must be a container of AlignedReads");
    const auto overlapped = overlap_range(reads, region);
    std::vector<double> qualities {};
    qualities.reserve(count_base_pairs(reads, region, NonMapTag {}));
    std::for_each(std::cbegin(overlapped), std::cend(overlapped), [&qualities] (const auto& read) {
        for (const auto quality : read.base_qualities()) {
            qualities.push_back(static_cast<double>(quality));
        }
    });
    return maths::median<double>(qualities);
}

template <typename T>
double rmq_base_quality(const T& reads, NonMapTag)
{
    static_assert(is_aligned_read_container<T>, "T must be a container of AlignedReads");
    std::vector<double> qualities {};
    qualities.reserve(count_base_pairs(reads, NonMapTag {}));
    for (const auto& read : reads) {
        for (const auto quality : read.base_qualities()) {
            qualities.push_back(static_cast<double>(quality));
        }
    }
    return maths::rmq<double>(qualities);
}

template <typename T>
double rmq_base_quality(const T& reads, const GenomicRegion& region, NonMapTag)
{
    static_assert(is_aligned_read_container<T>, "T must be a container of AlignedReads");
    const auto overlapped = overlap_range(reads, region);
    std::vector<double> qualities {};
    qualities.reserve(count_base_pairs(reads, region, NonMapTag {}));
    std::for_each(std::cbegin(overlapped), std::cend(overlapped), [&qualities] (const auto& read) {
        for (const auto quality : read.base_qualities()) {
            qualities.push_back(static_cast<double>(quality));
        }
    });
    return maths::rmq<double>(qualities);
}

template <typename ReadMap>
bool has_coverage(const ReadMap& reads, MapTag)
{
    return std::any_of(std::cbegin(reads), std::cend(reads),
                       [] (const auto& sample_reads) {
                           return has_coverage(sample_reads.second, NonMapTag {});
                       });
}

template <typename ReadMap>
bool has_coverage(const ReadMap& reads, const GenomicRegion& region, MapTag)
{
    return std::any_of(std::cbegin(reads), std::cend(reads),
                       [&region] (const auto& sample_reads) {
                           return has_coverage(sample_reads.second, region, NonMapTag {});
                       });
}

template <typename ReadMap>
unsigned min_coverage(const ReadMap& reads, MapTag)
{
    if (reads.empty()) return 0;
    std::vector<unsigned> sample_min_coverages(reads.size(), 0);
    std::transform(std::cbegin(reads), std::cend(reads), std::begin(sample_min_coverages),
                   [] (const auto& sample_reads) {
                       return min_coverage(sample_reads.second, NonMapTag {});
                   });
    return *std::min_element(cbegin(sample_min_coverages), std::cend(sample_min_coverages));
}

template <typename ReadMap>
unsigned min_coverage(const ReadMap& reads, const GenomicRegion& region, MapTag)
{
    if (reads.empty()) return 0;
    std::vector<unsigned> sample_min_coverages(reads.size(), 0);
    std::transform(std::cbegin(reads), std::cend(reads), std::begin(sample_min_coverages),
                   [&region] (const auto& sample_reads) {
                       return min_coverage(sample_reads.second, region, NonMapTag {});
                   });
    return *std::min_element(cbegin(sample_min_coverages), std::cend(sample_min_coverages));
}

template <typename ReadMap>
unsigned max_coverage(const ReadMap& reads, MapTag)
{
    if (reads.empty()) return 0;
    std::vector<unsigned> sample_max_coverages(reads.size(), 0);
    std::transform(std::cbegin(reads), std::cend(reads), std::begin(sample_max_coverages),
                   [] (const auto& sample_reads) {
                       return max_coverage(sample_reads.second, NonMapTag {});
                   });
    return *std::max_element(cbegin(sample_max_coverages), std::cend(sample_max_coverages));
}

template <typename ReadMap>
unsigned max_coverage(const ReadMap& reads, const GenomicRegion& region, MapTag)
{
    if (reads.empty()) return 0;
    std::vector<unsigned> sample_max_coverages(reads.size(), 0);
    std::transform(std::cbegin(reads), std::cend(reads), std::begin(sample_max_coverages),
                   [&region] (const auto& sample_reads) {
                       return max_coverage(sample_reads.second, region, NonMapTag {});
                   });
    return *std::max_element(std::cbegin(sample_max_coverages), std::cend(sample_max_coverages));
}

template <typename T>
auto sum_positional_coverages(const T& reads, const GenomicRegion& region)
{
    const auto num_bases = size(region);
    std::vector<std::size_t> result(num_bases);
    for (const auto& p : reads) {
        const auto sample_coverages = calculate_positional_coverage(p.second, region);
        assert(sample_coverages.size() == num_bases);
        for (std::size_t i {0}; i < num_bases; ++i) {
            result[i] += sample_coverages[i];
        }
    }
    return result;
}

template <typename T>
double mean_coverage(const T& reads, const GenomicRegion& region, MapTag)
{
    if (reads.empty() || is_empty_region(region)) return 0.0;
    if (reads.size() == 1) return mean_coverage(std::cbegin(reads)->second, region, NonMapTag {});
    return maths::mean(sum_positional_coverages(reads, region));
}

template <typename T>
double mean_coverage(const T& reads, MapTag)
{
    if (reads.empty() || all_empty(reads)) return 0.0;
    return mean_coverage(reads, encompassing_region(reads), MapTag {});
}

template <typename T>
double stdev_coverage(const T& reads, const GenomicRegion& region, MapTag)
{
    if (reads.empty() || is_empty_region(region)) return 0.0;
    if (reads.size() == 1) return stdev_coverage(std::cbegin(reads)->second, region, NonMapTag {});
    return maths::stdev(sum_positional_coverages(reads, region));
}

template <typename T>
double stdev_coverage(const T& reads, MapTag)
{
    if (reads.empty() || all_empty(reads)) return 0.0;
    return stdev_coverage(reads, encompassing_region(reads), MapTag {});
}

template <typename T>
std::size_t count_reads(const T& reads, MapTag)
{
    return maths::sum_sizes(reads);
}

template <typename T>
std::size_t count_reads(const T& reads, const GenomicRegion& region, MapTag)
{
    return std::accumulate(std::cbegin(reads), std::cend(reads), std::size_t {0},
                           [&region] (const auto curr, const auto& sample_reads) {
                               return curr + count_reads(sample_reads.second, region, NonMapTag {});
                           });
}

template <typename T>
std::size_t count_forward(const T& reads, MapTag)
{
    return std::accumulate(std::cbegin(reads), std::cend(reads), std::size_t {0},
                           [] (const auto curr, const auto& sample_reads) {
                               return curr + count_forward(sample_reads.second, NonMapTag {});
                           });
}

template <typename T>
std::size_t count_forward(const T& reads, const GenomicRegion& region, MapTag)
{
    return std::accumulate(std::cbegin(reads), std::cend(reads), std::size_t {0},
                           [&region] (const auto curr, const auto& sample_reads) {
                               return curr + count_forward(sample_reads.second, region, NonMapTag {});
                           });
}

template <typename T>
std::size_t count_reverse(const T& reads, MapTag)
{
    return std::accumulate(std::cbegin(reads), std::cend(reads), std::size_t {0},
                           [] (const auto curr, const auto& sample_reads) {
                               return curr + count_reverse(sample_reads.second, NonMapTag {});
                           });
}

template <typename T>
std::size_t count_reverse(const T& reads, const GenomicRegion& region, MapTag)
{
    return std::accumulate(std::cbegin(reads), std::cend(reads), std::size_t {0},
                           [&region] (const auto curr, const auto& sample_reads) {
                               return curr + count_reverse(sample_reads.second, region, NonMapTag {});
                           });
}

template <typename T>
std::pair<std::size_t, std::size_t> count_directions(const T& reads, MapTag)
{
    std::pair<std::size_t, std::size_t> result {};
    for (const auto& p : reads) {
        const auto sample_counts = count_directions(p.second, NonMapTag {});
        result.first += sample_counts.first; result.second += sample_counts.second;
    }
    return result;
}

template <typename T>
std::pair<std::size_t, std::size_t> count_directions(const T& reads, const GenomicRegion& region, MapTag)
{
    std::pair<std::size_t, std::size_t> result {};
    for (const auto& p : reads) {
        const auto sample_counts = count_directions(p.second, region, NonMapTag {});
        result.first += sample_counts.first; result.second += sample_counts.second;
    }
    return result;
}

template <typename T>
AlignedRead::NucleotideSequence::size_type min_read_length(const T& reads, MapTag)
{
    boost::optional<AlignedRead::NucleotideSequence::size_type> min_length {};
    for (const auto& p : reads) {
        if (!p.second.empty()) {
            auto l = min_read_length(p.second, NonMapTag {});
            if (min_length) {
                min_length = std::min(*min_length, l);
            } else {
                min_length = l;
            }
        }
    }
    if (min_length) {
        return *min_length;
    } else {
        return 0;
    }
}

template <typename T>
AlignedRead::NucleotideSequence::size_type min_read_length(const T& reads, const GenomicRegion& region, MapTag)
{
    boost::optional<AlignedRead::NucleotideSequence::size_type> min_length {};
    for (const auto& p : reads) {
        if (!p.second.empty()) {
            auto l = min_read_length(p.second, region, NonMapTag {});
            if (min_length) {
                min_length = std::min(*min_length, l);
            } else {
                min_length = l;
            }
        }
    }
    if (min_length) {
        return *min_length;
    } else {
        return 0;
    }
}

template <typename T>
AlignedRead::NucleotideSequence::size_type max_read_length(const T& reads, MapTag)
{
    boost::optional<AlignedRead::NucleotideSequence::size_type> max_length {};
    for (const auto& p : reads) {
        if (!p.second.empty()) {
            auto l = max_read_length(p.second, NonMapTag {});
            if (max_length) {
                max_length = std::min(*max_length, l);
            } else {
                max_length = l;
            }
        }
    }
    if (max_length) {
        return *max_length;
    } else {
        return 0;
    }
}

template <typename T>
AlignedRead::NucleotideSequence::size_type max_read_length(const T& reads, const GenomicRegion& region, MapTag)
{
    boost::optional<AlignedRead::NucleotideSequence::size_type> max_length {};
    for (const auto& p : reads) {
        if (!p.second.empty()) {
            auto l = max_read_length(p.second, region, NonMapTag {});
            if (max_length) {
                max_length = std::min(*max_length, l);
            } else {
                max_length = l;
            }
        }
    }
    if (max_length) {
        return *max_length;
    } else {
        return 0;
    }
}

template <typename T>
std::size_t count_base_pairs(const T& reads, MapTag)
{
    return std::accumulate(std::cbegin(reads), std::cend(reads), std::size_t {0},
                           [] (const auto curr, const auto& sample_reads) {
                               return curr + count_base_pairs(sample_reads.second, NonMapTag {});
                           });
}

template <typename T>
std::size_t count_base_pairs(const T& reads, const GenomicRegion& region, MapTag)
{
    return std::accumulate(std::cbegin(reads), std::cend(reads), std::size_t {0},
                           [&region] (const auto curr, const auto& sample_reads) {
                               return curr + count_base_pairs(sample_reads.second, region, NonMapTag {});
                           });
}

template <typename T>
std::size_t count_forward_base_pairs(const T& reads, MapTag)
{
    return std::accumulate(std::cbegin(reads), std::cend(reads), std::size_t {0},
                           [] (const auto curr, const auto& sample_reads) {
                               return curr + count_forward_base_pairs(sample_reads.second, NonMapTag {});
                           });
}

template <typename T>
std::size_t count_forward_base_pairs(const T& reads, const GenomicRegion& region, MapTag)
{
    return std::accumulate(std::cbegin(reads), std::cend(reads), std::size_t {0},
                           [&region] (const auto curr, const auto& sample_reads) {
                               return curr + count_forward_base_pairs(sample_reads.second, region, NonMapTag {});
                           });
}

template <typename T>
std::size_t count_reverse_base_pairs(const T& reads, MapTag)
{
    return std::accumulate(std::cbegin(reads), std::cend(reads), std::size_t {0},
                           [] (const auto curr, const auto& sample_reads) {
                               return curr + count_reverse_base_pairs(sample_reads.second, NonMapTag {});
                           });
}

template <typename T>
std::size_t count_reverse_base_pairs(const T& reads, const GenomicRegion& region, MapTag)
{
    return std::accumulate(std::cbegin(reads), std::cend(reads), std::size_t {0},
                           [&region] (const auto curr, const auto& sample_reads) {
                               return curr + count_reverse_base_pairs(sample_reads.second, region, NonMapTag {});
                           });
}

template <typename T>
std::size_t count_mapq_zero(const T& reads, MapTag)
{
    return std::accumulate(std::cbegin(reads), std::cend(reads), std::size_t {0},
                           [] (const auto curr, const auto& sample_reads) {
                               return curr + count_mapq_zero(sample_reads.second, NonMapTag {});
                           });
}

template <typename T>
std::size_t count_mapq_zero(const T& reads, const GenomicRegion& region, MapTag)
{
    return std::accumulate(std::cbegin(reads), std::cend(reads), std::size_t {0},
                           [&region] (const auto curr, const auto& sample_reads) {
                               return curr + count_mapq_zero(sample_reads.second, region, NonMapTag {});
                           });
}

template <typename T>
auto min_mapping_quality(const T& reads, MapTag)
{
    auto result = min_mapping_quality(std::cbegin(reads)->second, NonMapTag {});
    std::for_each(std::next(std::cbegin(reads)), std::cend(reads), [&] (const auto& p) {
        result = std::min(result, min_mapping_quality(p.second, NonMapTag {}));
    });
    return result;
}

template <typename T>
auto min_mapping_quality(const T& reads, const GenomicRegion& region, MapTag)
{
    auto result = min_mapping_quality(std::cbegin(reads)->second, region, NonMapTag {});
    std::for_each(std::next(std::cbegin(reads)), std::cend(reads), [&] (const auto& p) {
        result = std::min(result, min_mapping_quality(p.second, region, NonMapTag {}));
    });
    return result;
}

template <typename T>
auto max_mapping_quality(const T& reads, MapTag)
{
    auto result = max_mapping_quality(std::cbegin(reads)->second, NonMapTag {});
    std::for_each(std::next(std::cbegin(reads)), std::cend(reads), [&] (const auto& p) {
        result = std::max(result, max_mapping_quality(p.second, NonMapTag {}));
    });
    return result;
}

template <typename T>
auto max_mapping_quality(const T& reads, const GenomicRegion& region, MapTag)
{
    auto result = max_mapping_quality(std::cbegin(reads)->second, region, NonMapTag {});
    std::for_each(std::next(std::cbegin(reads)), std::cend(reads), [&] (const auto& p) {
        result = std::max(result, max_mapping_quality(p.second, region, NonMapTag {}));
    });
    return result;
}

template <typename T>
double median_mapping_quality(const T& reads, MapTag)
{
    std::vector<double> qualities {};
    qualities.reserve(maths::sum_sizes(reads));
    for (const auto& sample_reads : reads) {
        std::transform(std::cbegin(sample_reads.second), std::cend(sample_reads.second),
                       std::back_inserter(qualities), [] (const auto& read) {
            return static_cast<double>(read.mapping_quality());
        });
    }
    return maths::median<double>(qualities);
}

template <typename T>
double median_mapping_quality(const T& reads, const GenomicRegion& region, MapTag)
{
    std::vector<double> qualities {};
    qualities.reserve(maths::sum_sizes(reads));
    for (const auto& sample_reads : reads) {
        const auto overlapped = overlap_range(sample_reads.second, region);
        std::transform(std::cbegin(overlapped), std::cend(overlapped),
                       std::back_inserter(qualities), [] (const auto& read) {
            return static_cast<double>(read.mapping_quality());
        });
    }
    return maths::median<double>(qualities);
}

template <typename T>
double rmq_mapping_quality(const T& reads, MapTag)
{
    std::vector<double> qualities {};
    qualities.reserve(maths::sum_sizes(reads));
    for (const auto& sample_reads : reads) {
        std::transform(std::cbegin(sample_reads.second), std::cend(sample_reads.second),
                       std::back_inserter(qualities), [] (const auto& read) {
                           return static_cast<double>(read.mapping_quality());
                       });
    }
    return maths::rmq<double>(qualities);
}

template <typename T>
double rmq_mapping_quality(const T& reads, const GenomicRegion& region, MapTag)
{
    std::vector<double> qualities {};
    qualities.reserve(maths::sum_sizes(reads));
    for (const auto& sample_reads : reads) {
        const auto overlapped = overlap_range(sample_reads.second, region);
        std::transform(std::cbegin(overlapped), std::cend(overlapped),
                       std::back_inserter(qualities), [] (const auto& read) {
                           return static_cast<double>(read.mapping_quality());
                       });
    }
    return maths::rmq<double>(qualities);
}

template <typename T>
double median_base_quality(const T& reads, MapTag)
{
    std::vector<double> qualities {};
    qualities.reserve(count_base_pairs(reads, MapTag {}));
    for (const auto& sample_reads : reads) {
        for (const auto& read : sample_reads.second) {
            for (const auto quality : read.base_qualities()) {
                qualities.push_back(static_cast<double>(quality));
            }
        }
    }
    return maths::median<double>(qualities);
}

template <typename T>
double median_base_quality(const T& reads, const GenomicRegion& region, MapTag)
{
    std::vector<double> qualities {};
    qualities.reserve(count_base_pairs(reads, region, MapTag {}));
    for (const auto& sample_reads : reads) {
        const auto overlapped = overlap_range(sample_reads.second, region);
        std::for_each(std::cbegin(overlapped), std::cend(overlapped),
                      [&qualities] (const auto& read) {
                          for (const auto quality : read.base_qualities()) {
                              qualities.push_back(static_cast<double>(quality));
                          }
                      });
    }
    return maths::median<double>(qualities);
}

template <typename T>
double rmq_base_quality(const T& reads, MapTag)
{
    std::vector<double> qualities {};
    qualities.reserve(count_base_pairs(reads, MapTag {}));
    for (const auto& sample_reads : reads) {
        for (const auto& read : sample_reads.second) {
            for (const auto quality : read.base_qualities()) {
                qualities.push_back(static_cast<double>(quality));
            }
        }
    }
    return maths::rmq<double>(qualities);
}

template <typename T>
double rmq_base_quality(const T& reads, const GenomicRegion& region, MapTag)
{
    std::vector<double> qualities {};
    qualities.reserve(count_base_pairs(reads, region, MapTag {}));
    for (const auto& sample_reads : reads) {
        const auto overlapped = overlap_range(sample_reads.second, region);
        std::for_each(std::cbegin(overlapped), std::cend(overlapped),
                      [&qualities] (const auto& read) {
                          for (const auto quality : read.base_qualities()) {
                              qualities.push_back(static_cast<double>(quality));
                          }
                      });
    }
    return maths::rmq<double>(qualities);
}

} // namespace detail

template <typename T, typename = enable_if_map<T>>
bool has_coverage(const T& reads)
{
    return detail::has_coverage(reads, MapTag {});
}

template <typename T, typename = enable_if_map<T>>
bool has_coverage(const T& reads, const GenomicRegion& region)
{
    return detail::has_coverage(reads, region, MapTag {});
}

template <typename T, typename = enable_if_map<T>>
unsigned min_coverage(const T& reads)
{
    return detail::min_coverage(reads, MapTag {});
}

template <typename T, typename = enable_if_map<T>>
unsigned min_coverage(const T& reads, const GenomicRegion& region)
{
    return detail::min_coverage(reads, region, MapTag {});
}

template <typename T, typename = enable_if_map<T>>
unsigned max_coverage(const T& reads)
{
    return detail::max_coverage(reads, MapTag {});
}

template <typename T, typename = enable_if_map<T>>
unsigned max_coverage(const T& reads, const GenomicRegion& region)
{
    return detail::max_coverage(reads, region, MapTag {});
}

template <typename T>
double mean_coverage(const T& reads)
{
    return detail::mean_coverage(reads, MapTagType<T> {});
}

template <typename T>
double mean_coverage(const T& reads, const GenomicRegion& region)
{
    return detail::mean_coverage(reads, region, MapTagType<T> {});
}

template <typename T>
double stdev_coverage(const T& reads)
{
    return detail::stdev_coverage(reads, MapTagType<T> {});
}

template <typename T>
double stdev_coverage(const T& reads, const GenomicRegion& region)
{
    return detail::stdev_coverage(reads, region, MapTagType<T> {});
}

template <typename T>
std::size_t count_reads(const T& reads)
{
    return detail::count_reads(reads, MapTagType<T> {});
}

template <typename T>
std::size_t count_reads(const T& reads, const GenomicRegion& region)
{
    return detail::count_reads(reads, region, MapTagType<T> {});
}

template <typename T>
std::size_t count_forward(const T& reads)
{
    return detail::count_forward(reads, MapTagType<T> {});
}

template <typename T>
std::size_t count_forward(const T& reads, const GenomicRegion& region)
{
    return detail::count_forward(reads, region, MapTagType<T> {});
}

template <typename T>
std::size_t count_reverse(const T& reads)
{
    return detail::count_reverse(reads, MapTagType<T> {});
}

template <typename T>
std::size_t count_reverse(const T& reads, const GenomicRegion& region)
{
    return detail::count_reverse(reads, region, MapTagType<T> {});
}

template <typename T>
std::pair<std::size_t, std::size_t> count_directions(const T& reads)
{
    return detail::count_directions(reads, MapTagType<T> {});
}

template <typename T>
std::pair<std::size_t, std::size_t> count_directions(const T& reads, const GenomicRegion& region)
{
    return detail::count_directions(reads, region, MapTagType<T> {});
}

template <typename T>
AlignedRead::NucleotideSequence::size_type min_read_length(const T& reads)
{
    return detail::min_read_length(reads, MapTagType<T> {});
}

template <typename T>
AlignedRead::NucleotideSequence::size_type min_read_length(const T& reads, const GenomicRegion& region)
{
    return detail::min_read_length(reads, region, MapTagType<T> {});
}

template <typename T>
AlignedRead::NucleotideSequence::size_type max_read_length(const T& reads)
{
    return detail::max_read_length(reads, MapTagType<T> {});
}

template <typename T>
AlignedRead::NucleotideSequence::size_type max_read_length(const T& reads, const GenomicRegion& region)
{
    return detail::max_read_length(reads, region, MapTagType<T> {});
}

template <typename T>
double strand_bias(const T& reads)
{
    const auto num_forward_reads = static_cast<double>(count_forward(reads));
    const auto num_reverse_reads = static_cast<double>(count_reverse(reads));
    const auto total             = num_forward_reads + num_reverse_reads;
    return (total > 0) ? (num_forward_reads / total) : 0.0;
}

template <typename T>
double strand_bias(const T& reads, const GenomicRegion& region)
{
    const auto num_forward_reads = static_cast<double>(count_forward(reads, region));
    const auto num_reverse_reads = static_cast<double>(count_reverse(reads, region));
    const auto total             = num_forward_reads + num_reverse_reads;
    return (total > 0) ? (num_forward_reads / total) : 0.0;
}

template <typename T>
std::size_t count_base_pairs(const T& reads)
{
    return detail::count_base_pairs(reads, MapTagType<T> {});
}

template <typename T>
std::size_t count_base_pairs(const T& reads, const GenomicRegion& region)
{
    return detail::count_base_pairs(reads, region, MapTagType<T> {});
}

template <typename T>
std::size_t count_forward_base_pairs(const T& reads)
{
    return detail::count_forward_base_pairs(reads, MapTagType<T> {});
}

template <typename T>
std::size_t count_forward_base_pairs(const T& reads, const GenomicRegion& region)
{
    return detail::count_forward_base_pairs(reads, region, MapTagType<T> {});
}

template <typename T>
std::size_t count_reverse_base_pairs(const T& reads)
{
    return detail::count_reverse_base_pairs(reads, MapTagType<T> {});
}

template <typename T>
std::size_t count_reverse_base_pairs(const T& reads, const GenomicRegion& region)
{
    return detail::count_reverse_base_pairs(reads, region, MapTagType<T> {});
}

template <typename T>
std::size_t count_mapq_zero(const T& reads)
{
    return detail::count_mapq_zero(reads, MapTagType<T> {});
}

template <typename T>
std::size_t count_mapq_zero(const T& reads, const GenomicRegion& region)
{
    return detail::count_mapq_zero(reads, region, MapTagType<T> {});
}

template <typename T>
auto min_mapping_quality(const T& reads)
{
    return detail::min_mapping_quality(reads, MapTagType<T> {});
}

template <typename T>
auto min_mapping_quality(const T& reads, const GenomicRegion& region)
{
    return detail::min_mapping_quality(reads, region, MapTagType<T> {});
}

template <typename T>
auto max_mapping_quality(const T& reads)
{
    return detail::max_mapping_quality(reads, MapTagType<T> {});
}

template <typename T>
auto max_mapping_quality(const T& reads, const GenomicRegion& region)
{
    return detail::max_mapping_quality(reads, region, MapTagType<T> {});
}

template <typename T>
double median_mapping_quality(const T& reads)
{
    return detail::median_mapping_quality(reads, MapTagType<T> {});
}

template <typename T>
double median_mapping_quality(const T& reads, const GenomicRegion& region)
{
    return detail::median_mapping_quality(reads, region, MapTagType<T> {});
}

template <typename T>
double rmq_mapping_quality(const T& reads)
{
    return detail::rmq_mapping_quality(reads, MapTagType<T> {});
}

template <typename T>
double rmq_mapping_quality(const T& reads, const GenomicRegion& region)
{
    return detail::rmq_mapping_quality(reads, region, MapTagType<T> {});
}

template <typename T>
double median_base_quality(const T& reads)
{
    return detail::median_base_quality(reads, MapTagType<T> {});
}

template <typename T>
double median_base_quality(const T& reads, const GenomicRegion& region)
{
    return detail::median_base_quality(reads, region, MapTagType<T> {});
}

template <typename T>
double rmq_base_quality(const T& reads)
{
    return detail::rmq_base_quality(reads, MapTagType<T> {});
}

template <typename T>
double rmq_base_quality(const T& reads, const GenomicRegion& region)
{
    return detail::rmq_base_quality(reads, region, MapTagType<T> {});
}

template <typename ReadMap>
std::size_t count_samples_with_coverage(const ReadMap& reads)
{
    return std::count_if(std::cbegin(reads), std::cend(reads),
                         [] (const auto& sample_reads) {
                             return has_coverage(sample_reads.second);
                         });
}

template <typename ReadMap>
std::size_t count_samples_with_coverage(const ReadMap& reads, const GenomicRegion& region)
{
    return std::count_if(std::cbegin(reads), std::cend(reads),
                         [&region] (const auto& sample_reads) {
                             return has_coverage(sample_reads.second, region);
                         });
}

template <typename ReadMap>
unsigned sum_min_coverages(const ReadMap& reads)
{
    return std::accumulate(std::cbegin(reads), std::cend(reads), 0,
                           [] (auto curr, const auto& sample_reads) {
                               return curr + min_coverage(sample_reads.second);
                           });
}

template <typename ReadMap>
unsigned sum_min_coverages(const ReadMap& reads, const GenomicRegion& region)
{
    return std::accumulate(std::cbegin(reads), std::cend(reads), 0,
                           [&region] (auto curr, const auto& sample_reads) {
                               return curr + min_coverage(sample_reads.second, region);
                           });
}

template <typename ReadMap>
unsigned sum_max_coverages(const ReadMap& reads)
{
    return std::accumulate(std::cbegin(reads), std::cend(reads), 0,
                           [] (auto curr, const auto& sample_reads) {
                               return curr + max_coverage(sample_reads.second);
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

template <typename ReadMap>
std::size_t max_sample_read_count(const ReadMap& reads)
{
    if (reads.empty()) return 0;
    return std::max_element(std::cbegin(reads), std::cend(reads),
                            [] (const auto& lhs, const auto& rhs) {
                                return lhs.second.size() < rhs.second.size();
                            })->second.size();
}

} // namespace octopus

#endif
