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
#include <iterator>
#include <utility>
#include <algorithm>
#include <numeric>
#include <type_traits>

#include <boost/range/algorithm.hpp>

#include "aligned_read.hpp"
#include "mappable_algorithms.hpp"
#include "mappable_set.hpp"
#include "mappable_map.hpp"
#include "maths.hpp"
#include "type_tricks.hpp"

namespace Octopus
{

namespace detail
{
    namespace
    {
        template <typename T>
        constexpr bool is_aligned_read = std::is_same<std::decay_t<T>, AlignedRead>::value;
        
        template <typename Container>
        constexpr bool is_aligned_read_container = is_aligned_read<typename Container::value_type>;
        
        struct IsForward
        {
            bool operator()(const AlignedRead& read) const
            {
                return !read.is_marked_reverse_mapped();
            }
        };
        
        struct IsReverse
        {
            bool operator()(const AlignedRead& read) const
            {
                return read.is_marked_reverse_mapped();
            }
        };
        
        struct IsMappingQualityZero
        {
            bool operator()(const AlignedRead& read) const noexcept
            {
                return read.get_mapping_quality() == 0;
            }
        };
    } // namespace
    
    template <typename T>
    bool has_coverage(const T& reads, const GenomicRegion& region, NonMapTag)
    {
        const auto overlapped = overlap_range(reads, region);
        return std::any_of(std::cbegin(overlapped), std::cend(overlapped),
                           [] (const auto& read) { return !is_empty_region(read); });
    }
    
    inline bool has_coverage(const MappableSet<AlignedRead>& reads, NonMapTag)
    {
        return std::any_of(std::cbegin(reads), std::cend(reads),
                           [] (const auto& read) { return !is_empty_region(read); });
    }
    
    template <typename T>
    unsigned min_coverage(const T& reads, const GenomicRegion& region, NonMapTag)
    {
        if (reads.empty() || is_empty_region(region)) return 0;
        const auto positional_coverage = calculate_positional_coverage(reads, region);
        return *std::min_element(std::cbegin(positional_coverage), std::cend(positional_coverage));
    }
    
    inline unsigned min_coverage(const MappableSet<AlignedRead>& reads, NonMapTag)
    {
        if (reads.empty()) return 0;
        const auto positional_coverage = calculate_positional_coverage(reads);
        return *std::min_element(std::cbegin(positional_coverage), std::cend(positional_coverage));
    }
    
    template <typename T>
    unsigned max_coverage(const T& reads, const GenomicRegion& region, NonMapTag)
    {
        if (reads.empty() || is_empty_region(region)) return 0;
        const auto positional_coverage = calculate_positional_coverage(reads, region);
        return *std::max_element(std::cbegin(positional_coverage), std::cend(positional_coverage));
    }
    
    inline unsigned max_coverage(const MappableSet<AlignedRead>& reads, NonMapTag)
    {
        if (reads.empty()) return 0;
        const auto positional_coverage = calculate_positional_coverage(reads);
        return *std::max_element(std::cbegin(positional_coverage), std::cend(positional_coverage));
    }
    
    template <typename T>
    double mean_coverage(const T& reads, const GenomicRegion& region, NonMapTag)
    {
        if (reads.empty() || is_empty_region(region)) return 0;
        return Maths::mean(positional_coverage(reads, region));
    }
    
    inline double mean_coverage(const MappableSet<AlignedRead>& reads, NonMapTag)
    {
        if (reads.empty()) return 0;
        return Maths::mean(calculate_positional_coverage(reads));
    }
    
    template <typename T>
    double stdev_coverage(const T& reads, const GenomicRegion& region, NonMapTag)
    {
        if (reads.empty() || is_empty_region(region)) return 0;
        return Maths::stdev(positional_coverage(reads, region));
    }
    
    inline double stdev_coverage(const MappableSet<AlignedRead>& reads, NonMapTag)
    {
        if (reads.empty()) return 0;
        return Maths::stdev(calculate_positional_coverage(reads));
    }
    
    template <typename T>
    size_t count_reads(const T& reads, NonMapTag)
    {
        return reads.size();
    }
    
    template <typename T>
    size_t count_reads(const T& reads, const GenomicRegion& region, NonMapTag)
    {
        return count_overlapped(reads, region);
    }
    
    inline size_t count_reads(const MappableSet<AlignedRead>& reads, const GenomicRegion& region, NonMapTag)
    {
        return reads.count_overlapped(region);
    }
    
    template <typename T>
    size_t count_forward(const T& reads, NonMapTag)
    {
        static_assert(is_aligned_read_container<T>, "T must be a container of AlignedReads");
        return std::count_if(std::cbegin(reads), std::cend(reads), IsForward {});
    }
    
    template <typename T>
    size_t count_forward(const T& reads, const GenomicRegion& region, NonMapTag)
    {
        static_assert(is_aligned_read_container<T>, "T must be a container of AlignedReads");
        const auto overlapped = overlap_range(reads, region);
        return std::count_if(std::cbegin(overlapped), std::cend(overlapped), IsForward {});
    }
    
    inline size_t count_forward(const MappableSet<AlignedRead>& reads, const GenomicRegion& region, NonMapTag)
    {
        const auto overlapped = reads.overlap_range(region);
        return std::count_if(std::cbegin(overlapped), std::cend(overlapped), IsForward {});
    }
    
    template <typename T>
    size_t count_reverse(const T& reads, NonMapTag)
    {
        static_assert(is_aligned_read_container<T>, "T must be a container of AlignedReads");
        return std::count_if(std::cbegin(reads), std::cend(reads), IsReverse {});
    }
    
    template <typename T>
    size_t count_reverse(const T& reads, const GenomicRegion& region, NonMapTag)
    {
        static_assert(is_aligned_read_container<T>, "T must be a container of AlignedReads");
        const auto overlapped = overlap_range(reads, region);
        return std::count_if(std::cbegin(overlapped), std::cend(overlapped), IsReverse {});
    }
    
    inline size_t count_reverse(const MappableSet<AlignedRead>& reads, const GenomicRegion& region, NonMapTag)
    {
        const auto overlapped = reads.overlap_range(region);
        return std::count_if(std::cbegin(overlapped), std::cend(overlapped), IsReverse {});
    }
    
    template <typename T>
    size_t count_base_pairs(const T& reads, NonMapTag)
    {
        static_assert(is_aligned_read_container<T>, "T must be a container of AlignedReads");
        return std::accumulate(std::cbegin(reads), std::cend(reads), size_t {},
                               [] (const auto curr, const auto& read) {
                                   return curr + static_cast<size_t>(read.get_sequence_size());
                               });
    }
    
    template <typename T>
    size_t count_base_pairs(const T& reads, const GenomicRegion& region, NonMapTag)
    {
        static_assert(is_aligned_read_container<T>, "T must be a container of AlignedReads");
        const auto overlapped = overlap_range(reads, region);
        return std::accumulate(std::cbegin(overlapped), std::cend(overlapped), size_t {},
                               [&region] (const auto curr, const auto& read) {
                                   return curr + count_overlapped_bases(read, region);
                               });
    }
    
    template <typename T>
    size_t count_forward_base_pairs(const T& reads, NonMapTag)
    {
        static_assert(is_aligned_read_container<T>, "T must be a container of AlignedReads");
        return std::accumulate(std::cbegin(reads), std::cend(reads), size_t {},
                               [] (const auto curr, const auto& read) {
                                   return curr + ((IsForward()(read)) ? static_cast<size_t>(read.get_sequence_size()) : 0);
                               });
    }
    
    template <typename T>
    size_t count_forward_base_pairs(const T& reads, const GenomicRegion& region, NonMapTag)
    {
        static_assert(is_aligned_read_container<T>, "T must be a container of AlignedReads");
        const auto overlapped = overlap_range(reads, region);
        return std::accumulate(std::cbegin(overlapped), std::cend(overlapped), size_t {},
                               [&region] (const auto curr, const auto& read) {
                                   return curr + ((IsForward()(read)) ? static_cast<size_t>(num_overlapped_bases(read, region)) : 0);
                               });
    }
    
    template <typename T>
    size_t count_reverse_base_pairs(const T& reads, NonMapTag)
    {
        static_assert(is_aligned_read_container<T>, "T must be a container of AlignedReads");
        return std::accumulate(std::cbegin(reads), std::cend(reads), size_t {},
                               [] (const auto curr, const auto& read) {
                                   return curr + ((IsReverse()(read)) ? static_cast<size_t>(read.get_sequence_size()) : 0);
                               });
    }
    
    template <typename T>
    size_t count_reverse_base_pairs(const T& reads, const GenomicRegion& region, NonMapTag)
    {
        static_assert(is_aligned_read_container<T>, "T must be a container of AlignedReads");
        const auto overlapped = overlap_range(reads, region);
        return std::accumulate(std::cbegin(overlapped), std::cend(overlapped), size_t {},
                               [&region] (const auto curr, const auto& read) {
                                   return curr + ((IsReverse()(read)) ? static_cast<size_t>(num_overlapped_bases(read, region)) : 0);
                               });
    }
    
    template <typename T>
    size_t count_mapq_zero(const T& reads, NonMapTag)
    {
        static_assert(is_aligned_read_container<T>, "T must be a container of AlignedReads");
        return std::count_if(std::cbegin(reads), std::cend(reads), IsMappingQualityZero {});
    }
    
    template <typename T>
    size_t count_mapq_zero(const T& reads, const GenomicRegion& region, NonMapTag)
    {
        static_assert(is_aligned_read_container<T>, "T must be a container of AlignedReads");
        const auto overlapped = overlap_range(reads, region);
        return std::count_if(std::cbegin(overlapped), std::cend(overlapped), IsMappingQualityZero {});
    }
    
    inline size_t count_mapq_zero(const MappableSet<AlignedRead>& reads, const GenomicRegion& region, NonMapTag)
    {
        const auto overlapped = reads.overlap_range(region);
        return std::count_if(std::cbegin(overlapped), std::cend(overlapped), IsMappingQualityZero {});
    }
    
    template <typename T>
    double rmq_mapping_quality(const T& reads, NonMapTag)
    {
        static_assert(is_aligned_read_container<T>, "T must be a container of AlignedReads");
        
        std::vector<double> qualities(reads.size());
        
        std::transform(std::cbegin(reads), std::cend(reads), std::begin(qualities),
                       [] (const auto& read) { return static_cast<double>(read.get_mapping_quality()); });
        
        return Maths::rmq<double>(qualities);
    }
    
    template <typename T>
    double rmq_mapping_quality(const T& reads, const GenomicRegion& region, NonMapTag)
    {
        static_assert(is_aligned_read_container<T>, "T must be a container of AlignedReads");
        
        const auto overlapped = overlap_range(reads, region);
        
        std::vector<double> qualities(size(overlapped));
        
        std::transform(std::cbegin(overlapped), std::cend(overlapped), std::begin(qualities),
                       [] (const auto& read) { return static_cast<double>(read.get_mapping_quality()); });
        
        return Maths::rmq<double>(qualities);
    }
    
    inline double rmq_mapping_quality(const MappableSet<AlignedRead>& reads, const GenomicRegion& region, NonMapTag)
    {
        const auto overlapped = reads.overlap_range(region);
        
        std::vector<double> qualities(size(overlapped));
        
        std::transform(std::cbegin(overlapped), std::cend(overlapped), std::begin(qualities),
                       [] (const auto& read) { return static_cast<double>(read.get_mapping_quality()); });
        
        return Maths::rmq<double>(qualities);
    }
    
    template <typename T>
    double rmq_base_quality(const T& reads, NonMapTag)
    {
        static_assert(is_aligned_read_container<T>, "T must be a container of AlignedReads");
        
        std::vector<double> qualities {};
        qualities.reserve(count_base_pairs(reads, NonMapTag {}));
        
        for (const auto& read : reads) {
            for (const auto quality : read.get_qualities()) {
                qualities.push_back(static_cast<double>(quality));
            }
        }
        
        return Maths::rmq<double>(qualities);
    }
    
    template <typename T>
    double rmq_base_quality(const T& reads, const GenomicRegion& region, NonMapTag)
    {
        static_assert(is_aligned_read_container<T>, "T must be a container of AlignedReads");
        
        const auto overlapped = overlap_range(reads, region);
        
        std::vector<double> qualities {};
        qualities.reserve(count_base_pairs(reads, region, NonMapTag {}));
        
        std::for_each(std::cbegin(overlapped), std::cend(overlapped), [&qualities] (const auto& read) {
            for (const auto quality : read.get_qualities()) {
                qualities.push_back(static_cast<double>(quality));
            }
        });
        
        return Maths::rmq<double>(qualities);
    }
    
    inline double rmq_base_quality(const MappableSet<AlignedRead>& reads, const GenomicRegion& region, NonMapTag)
    {
        const auto overlapped = reads.overlap_range(region);
        
        std::vector<double> qualities {};
        qualities.reserve(count_base_pairs(reads, region, NonMapTag {}));
        
        std::for_each(std::cbegin(overlapped), std::cend(overlapped), [&qualities] (const auto& read) {
            for (const auto quality : read.get_qualities()) {
                qualities.push_back(static_cast<double>(quality));
            }
        });
        
        return Maths::rmq<double>(qualities);
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
    double mean_coverage(const T& reads, MapTag)
    {
        if (reads.empty()) return 0.0;
        
        std::vector<double> sample_mean_coverages(reads.size(), 0.0);
        
        std::transform(std::cbegin(reads), std::cend(reads), std::begin(sample_mean_coverages),
                       [] (const auto& sample_reads) {
                           return mean_coverage(sample_reads.second, NonMapTag {});
                       });
        
        return Maths::mean(sample_mean_coverages);
    }
    
    template <typename T>
    double mean_coverage(const T& reads, const GenomicRegion& region, MapTag)
    {
        if (reads.empty()) return 0.0;
        
        std::vector<double> sample_mean_coverages(reads.size(), 0.0);
        
        std::transform(std::cbegin(reads), std::cend(reads), std::begin(sample_mean_coverages),
                       [&region] (const auto& sample_reads) {
                           return mean_coverage(sample_reads.second, region, NonMapTag {});
                       });
        
        return Maths::mean(sample_mean_coverages);
    }
    
    template <typename T>
    double stdev_coverage(const T& reads, MapTag)
    {
        if (reads.empty()) return 0.0;
        
        std::vector<double> sample_stdev_coverages(reads.size(), 0.0);
        
        std::transform(std::cbegin(reads), std::cend(reads), std::begin(sample_stdev_coverages),
                       [] (const auto& sample_reads) {
                           return stdev_coverage(sample_reads.second, NonMapTag {});
                       });
        
        return Maths::stdev(sample_stdev_coverages);
    }
    
    template <typename T>
    double stdev_coverage(const T& reads, const GenomicRegion& region, MapTag)
    {
        if (reads.empty()) return 0.0;
        
        std::vector<double> sample_stdev_coverages(reads.size(), 0.0);
        
        std::transform(std::cbegin(reads), std::cend(reads), std::begin(sample_stdev_coverages),
                       [&region] (const auto& sample_reads) {
                           return stdev_coverage(sample_reads.second, region, NonMapTag {});
                       });
        
        return Maths::stdev(sample_stdev_coverages);
    }
    
    template <typename T>
    size_t count_reads(const T& reads, MapTag)
    {
        return Maths::sum_sizes(reads);
    }
    
    template <typename T>
    size_t count_reads(const T& reads, const GenomicRegion& region, MapTag)
    {
        return std::accumulate(std::cbegin(reads), std::cend(reads), size_t {},
                               [&region] (const auto curr, const auto& sample_reads) {
                                   return curr + count_reads(sample_reads.second, region, NonMapTag {});
                               });
    }
    
    template <typename T>
    size_t count_forward(const T& reads, MapTag)
    {
        return std::accumulate(std::cbegin(reads), std::cend(reads), size_t {},
                               [] (const auto curr, const auto& sample_reads) {
                                   return curr + count_forward(sample_reads.second, NonMapTag {});
                               });
    }
    
    template <typename T>
    size_t count_forward(const T& reads, const GenomicRegion& region, MapTag)
    {
        return std::accumulate(std::cbegin(reads), std::cend(reads), size_t {},
                               [&region] (const auto curr, const auto& sample_reads) {
                                   return curr + count_forward(sample_reads.second, region, NonMapTag {});
                               });
    }
    
    template <typename T>
    size_t count_reverse(const T& reads, MapTag)
    {
        return std::accumulate(std::cbegin(reads), std::cend(reads), size_t {},
                               [] (const auto curr, const auto& sample_reads) {
                                   return curr + count_reverse(sample_reads.second, NonMapTag {});
                               });
    }
    
    template <typename T>
    size_t count_reverse(const T& reads, const GenomicRegion& region, MapTag)
    {
        return std::accumulate(std::cbegin(reads), std::cend(reads), size_t {},
                               [&region] (const auto curr, const auto& sample_reads) {
                                   return curr + count_reverse(sample_reads.second, region, NonMapTag {});
                               });
    }
    
    template <typename T>
    size_t count_base_pairs(const T& reads, MapTag)
    {
        return std::accumulate(std::cbegin(reads), std::cend(reads), size_t {},
                               [] (const auto curr, const auto& sample_reads) {
                                   return curr + count_base_pairs(sample_reads.second, NonMapTag {});
                               });
    }
    
    template <typename T>
    size_t count_base_pairs(const T& reads, const GenomicRegion& region, MapTag)
    {
        return std::accumulate(std::cbegin(reads), std::cend(reads), size_t {},
                               [&region] (const auto curr, const auto& sample_reads) {
                                   return curr + count_base_pairs(sample_reads.second, region, NonMapTag {});
                               });
    }
    
    template <typename T>
    size_t count_forward_base_pairs(const T& reads, MapTag)
    {
        return std::accumulate(std::cbegin(reads), std::cend(reads), size_t {},
                               [] (const auto curr, const auto& sample_reads) {
                                   return curr + count_forward_base_pairs(sample_reads.second, NonMapTag {});
                               });
    }
    
    template <typename T>
    size_t count_forward_base_pairs(const T& reads, const GenomicRegion& region, MapTag)
    {
        return std::accumulate(std::cbegin(reads), std::cend(reads), size_t {},
                               [&region] (const auto curr, const auto& sample_reads) {
                                   return curr + count_forward_base_pairs(sample_reads.second, region, NonMapTag {});
                               });
    }
    
    template <typename T>
    size_t count_reverse_base_pairs(const T& reads, MapTag)
    {
        return std::accumulate(std::cbegin(reads), std::cend(reads), size_t {},
                               [] (const auto curr, const auto& sample_reads) {
                                   return curr + count_reverse_base_pairs(sample_reads.second, NonMapTag {});
                               });
    }
    
    template <typename T>
    size_t count_reverse_base_pairs(const T& reads, const GenomicRegion& region, MapTag)
    {
        return std::accumulate(std::cbegin(reads), std::cend(reads), size_t {},
                               [&region] (const auto curr, const auto& sample_reads) {
                                   return curr + count_reverse_base_pairs(sample_reads.second, region, NonMapTag {});
                               });
    }
    
    template <typename T>
    size_t count_mapq_zero(const T& reads, MapTag)
    {
        return std::accumulate(std::cbegin(reads), std::cend(reads), size_t {},
                               [] (const auto curr, const auto& sample_reads) {
                                   return curr + count_mapq_zero(sample_reads.second, NonMapTag {});
                               });
    }
    
    template <typename T>
    size_t count_mapq_zero(const T& reads, const GenomicRegion& region, MapTag)
    {
        return std::accumulate(std::cbegin(reads), std::cend(reads), size_t {},
                               [&region] (const auto curr, const auto& sample_reads) {
                                   return curr + count_mapq_zero(sample_reads.second, region, NonMapTag {});
                               });
    }
    
    template <typename T>
    double rmq_mapping_quality(const T& reads, MapTag)
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
    double rmq_mapping_quality(const T& reads, const GenomicRegion& region, MapTag)
    {
        std::vector<double> qualities {};
        qualities.reserve(Maths::sum_sizes(reads));
        
        for (const auto& sample_reads : reads) {
            const auto overlapped = overlap_range(sample_reads.second, region);
            
            std::transform(std::cbegin(overlapped), std::cend(overlapped),
                           std::back_inserter(qualities), [] (const auto& read) {
                               return static_cast<double>(read.get_mapping_quality());
                           });
        }
        
        return Maths::rmq<double>(qualities);
    }
    
    template <typename T>
    double rmq_base_quality(const T& reads, MapTag)
    {
        std::vector<double> qualities {};
        qualities.reserve(count_base_pairs(reads, MapTag {}));
        
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
    double rmq_base_quality(const T& reads, const GenomicRegion& region, MapTag)
    {
        std::vector<double> qualities {};
        qualities.reserve(count_base_pairs(reads, region, MapTag {}));
        
        for (const auto& sample_reads : reads) {
            const auto overlapped = overlap_range(sample_reads.second, region);
            
            std::for_each(std::cbegin(overlapped), std::cend(overlapped),
                          [&qualities] (const auto& read) {
                              for (const auto quality : read.get_qualities()) {
                                  qualities.push_back(static_cast<double>(quality));
                              }
                          });
        }
        
        return Maths::rmq<double>(qualities);
    }
} // namespace detail

template <typename T>
bool has_coverage(const T& reads)
{
    return detail::has_coverage(reads, MapTagType<T> {});
}

template <typename T>
bool has_coverage(const T& reads, const GenomicRegion& region)
{
    return detail::has_coverage(reads, region, MapTagType<T> {});
}

template <typename T>
unsigned min_coverage(const T& reads)
{
    return detail::min_coverage(reads, MapTagType<T> {});
}

template <typename T>
unsigned min_coverage(const T& reads, const GenomicRegion& region)
{
    return detail::min_coverage(reads, region, MapTagType<T> {});
}

template <typename T>
unsigned max_coverage(const T& reads)
{
    return detail::max_coverage(reads, MapTagType<T> {});
}

template <typename T>
unsigned max_coverage(const T& reads, const GenomicRegion& region)
{
    return detail::max_coverage(reads, region, MapTagType<T> {});
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
size_t count_reads(const T& reads)
{
    return detail::count_reads(reads, MapTagType<T> {});
}

template <typename T>
size_t count_reads(const T& reads, const GenomicRegion& region)
{
    return detail::count_reads(reads, region, MapTagType<T> {});
}

template <typename T>
size_t count_forward(const T& reads)
{
    return detail::count_forward(reads, MapTagType<T> {});
}

template <typename T>
size_t count_forward(const T& reads, const GenomicRegion& region)
{
    return detail::count_forward(reads, region, MapTagType<T> {});
}

template <typename T>
size_t count_reverse(const T& reads)
{
    return detail::count_reverse(reads, MapTagType<T> {});
}

template <typename T>
size_t count_reverse(const T& reads, const GenomicRegion& region)
{
    return detail::count_reverse(reads, region, MapTagType<T> {});
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
size_t count_base_pairs(const T& reads)
{
    return detail::count_base_pairs(reads, MapTagType<T> {});
}

template <typename T>
size_t count_base_pairs(const T& reads, const GenomicRegion& region)
{
    return detail::count_base_pairs(reads, region, MapTagType<T> {});
}

template <typename T>
size_t count_forward_base_pairs(const T& reads)
{
    return detail::count_forward_base_pairs(reads, MapTagType<T> {});
}

template <typename T>
size_t count_forward_base_pairs(const T& reads, const GenomicRegion& region)
{
    return detail::count_forward_base_pairs(reads, region, MapTagType<T> {});
}

template <typename T>
size_t count_reverse_base_pairs(const T& reads)
{
    return detail::count_reverse_base_pairs(reads, MapTagType<T> {});
}

template <typename T>
size_t count_reverse_base_pairs(const T& reads, const GenomicRegion& region)
{
    return detail::count_reverse_base_pairs(reads, region, MapTagType<T> {});
}

template <typename T>
size_t count_mapq_zero(const T& reads)
{
    return detail::count_mapq_zero(reads, MapTagType<T> {});
}

template <typename T>
size_t count_mapq_zero(const T& reads, const GenomicRegion& region)
{
    return detail::count_mapq_zero(reads, region, MapTagType<T> {});
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
size_t count_samples_with_coverage(const ReadMap& reads)
{
    return std::count_if(std::cbegin(reads), std::cend(reads),
                         [] (const auto& sample_reads) {
                             return has_coverage(sample_reads.second);
                         });
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
size_t max_sample_read_count(const ReadMap& reads)
{
    if (reads.empty()) return 0;
    return std::max_element(std::cbegin(reads), std::cend(reads),
                            [] (const auto& lhs, const auto& rhs) {
                                return lhs.second.size() < rhs.second.size();
                            })->second.size();
}

namespace detail
{
    template <typename T, typename F>
    std::unordered_map<AlignedRead, unsigned>
    coverages_in_read_regions(const T& reads, const GenomicRegion& region, F f) {
        std::unordered_map<AlignedRead, unsigned> result {};
        result.reserve(reads.size());
        
        const auto position_coverages = positional_coverage(reads, region);
        
        const auto first_position = region_begin(region);
        const auto num_positions  = region_size(region);
        
        for (const auto& read : reads) {
            auto first = std::next(std::cbegin(position_coverages), (region_begin(read) <= first_position) ? 0 : region_begin(read) - first_position);
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
    
    auto positional_coverage = calculate_positional_coverage(reads, region);
    
    using Iterator = typename decltype(positional_coverage)::const_iterator;
    
    Iterator first {positional_coverage.cbegin()};
    Iterator current {first};
    Iterator last {positional_coverage.cend()};
    Iterator high_range_first, high_range_last;
    
    SizeType high_range_begin, high_range_end;
    
    while (current != last) {
        auto is_high_coverage = [max_coverage] (unsigned coverage) { return coverage > max_coverage; };
        
        high_range_first = std::find_if(current, last, is_high_coverage);
        
        if (high_range_first == last) break;
        
        high_range_last = std::find_if_not(high_range_first, last, is_high_coverage);
        
        high_range_begin = region_begin(region) + static_cast<SizeType>(std::distance(first, high_range_first));
        high_range_end   = high_range_begin + static_cast<SizeType>(std::distance(high_range_first, high_range_last));
        
        result.emplace_back(contig_name(region), high_range_begin, high_range_end);
        
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

template <typename T>
std::vector<GenomicRegion> find_uniform_coverage_regions(const T& reads, const GenomicRegion& region)
{
    const auto coverages = positional_coverage(reads, region);
    
    std::vector<GenomicRegion> result {};
    result.reserve(region_size(region));
    
    if (coverages.empty()) return result;
    
    const auto contig = contig_name(region);
    
    auto begin = region_begin(region);
    auto end   = begin;
    
    auto previous_coverage = coverages.front();
    
    for (const auto coverage : coverages) {
        if (coverage != previous_coverage) {
            result.emplace_back(contig, begin, end);
            begin = end;
            previous_coverage = coverage;
        }
        ++end;
    }
    
    result.emplace_back(contig, begin, end);
    
    result.shrink_to_fit();
    
    return result;
}

template <typename T>
std::vector<GenomicRegion> find_uniform_coverage_regions(const T& reads)
{
    return find_uniform_coverage_regions(reads, get_encompassing_region(reads));
}
    
namespace detail {
    inline MappableSet<AlignedRead> splice_all(const MappableSet<AlignedRead>& reads, const GenomicRegion& region, NonMapTag)
    {
        MappableSet<AlignedRead> result {};
        result.reserve(reads.size());
        
        for (const auto& read : reads) {
            result.emplace(splice(read, region));
        }
        
        return result;
    }
    
    template <typename T>
    T splice_all(const T& reads, const GenomicRegion& region, NonMapTag)
    {
        T result {};
        result.reserve(reads.size());
        
        std::transform(std::cbegin(reads), std::cend(reads), std::back_inserter(result),
                       [&region] (const auto& read) { return splice(read, region); });
        
        return result;
    }
    
    template <typename T>
    T splice_all(const T& reads, const GenomicRegion& region, MapTag)
    {
        T result {};
        result.reserve(reads.size());
        
        for (const auto& p : reads) {
            result.emplace(p.first, splice_all(p.second, region, NonMapTag {}));
        }
        
        return result;
    }
} // namespace detail
    
template <typename T>
T splice_all(const T& reads, const GenomicRegion& region)
{
    return detail::splice_all(reads, region, MapTagType<T> {});
}

template <typename Reads>
void compress_reads(Reads& reads)
{
    for (auto& read : reads) {
        read.compress();
    }
}

template <typename Reads>
void decompress_reads(Reads& reads)
{
    for (auto& read : reads) {
        read.decompress();
    }
}

// TODO
AlignedRead find_next_segment(const AlignedRead& read, const MappableMap<GenomicRegion::ContigNameType, AlignedRead>& reads);

// TODO
MappableSet<AlignedRead> find_chimeras(const AlignedRead& read, const MappableSet<AlignedRead>& reads);

} // namespace Octopus

#endif /* defined(__Octopus__read_utils__) */
