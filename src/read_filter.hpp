//
//  read_filter.hpp
//  Octopus
//
//  Created by Daniel Cooke on 06/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__read_filter__
#define __Octopus__read_filter__

#include <vector>
#include <memory>
#include <algorithm>
#include <iterator>
#include <type_traits>
#include <utility>
#include <unordered_map>

#include <boost/optional.hpp>

#include "common.hpp"
#include "read_filters.hpp"
#include "mappable_flat_multi_set.hpp"
#include "mappable_map.hpp"
#include "logging.hpp"

namespace Octopus
{
/*
 ReadFilter stores a collection of filter functions, which can either be non-context-based (
 applied to a single read), or context-based (applied to a range of reads).
 
 This is a template class as the type of iterator used for context-based filteration needs to be 
 known at compile time. The class needs to know what container it is going to be operating on.
 */
template <typename BidirIt_>
class ReadFilter
{
public:
    using BidirIt = BidirIt_;
    
    using BasicFilterPtr   = std::unique_ptr<ReadFilters::BasicReadFilter>;
    using ContextFilterPtr = std::unique_ptr<ReadFilters::ContextReadFilter<BidirIt>>;
    
    using FilterCountMap = std::unordered_map<std::string, std::size_t>;
    
    ReadFilter()  = default;
    ~ReadFilter() = default;
    
    ReadFilter(const ReadFilter&)            = default;
    ReadFilter& operator=(const ReadFilter&) = default;
    ReadFilter(ReadFilter&&)                 = default;
    ReadFilter& operator=(ReadFilter&&)      = default;
    
    void register_filter(BasicFilterPtr filter);
    void register_filter(ContextFilterPtr filter);
    
    unsigned num_filters() const noexcept;
    
    // These methods work just like std::remove and std::partition
    
    BidirIt remove(BidirIt first, BidirIt last) const;
    BidirIt partition(BidirIt first, BidirIt last) const;
    
    BidirIt remove(BidirIt first, BidirIt last, FilterCountMap& filter_counts) const;
    BidirIt partition(BidirIt first, BidirIt last, FilterCountMap& filter_counts) const;
    
private:
    std::vector<BasicFilterPtr> basic_filters_;
    std::vector<ContextFilterPtr> context_filters_;
    
    bool passes_all_basic_filters(const AlignedRead& read) const noexcept;
    auto find_failing_basic_filter(const AlignedRead& read) const noexcept;
};

template <typename BidirIt>
void ReadFilter<BidirIt>::register_filter(BasicFilterPtr filter)
{
    basic_filters_.emplace_back(std::move(filter));
}

template <typename BidirIt>
void ReadFilter<BidirIt>::register_filter(ContextFilterPtr filter)
{
    context_filters_.emplace_back(std::move(filter));
}

template <typename BidirIt>
unsigned ReadFilter<BidirIt>::num_filters() const noexcept
{
    return static_cast<unsigned>(basic_filters_.size() + context_filters_.size());
}

template <typename BidirIt>
BidirIt ReadFilter<BidirIt>::remove(BidirIt first, BidirIt last) const
{
    if (first == last || num_filters() == 0) return last;
    
    if (!basic_filters_.empty()) {
        last = std::remove_if(first, last,
                              [this] (const auto& read) {
                                  return !passes_all_basic_filters(read);
                              });
    }
    
    std::for_each(cbegin(context_filters_), cend(context_filters_),
                  [first, &last] (const auto& filter) {
                      last = filter->remove(first, last);
                  });
    
    return last;
}

template <typename BidirIt>
BidirIt ReadFilter<BidirIt>::partition(BidirIt first, BidirIt last) const
{
    if (first == last || num_filters() == 0) return last;
    
    const auto passes_basic_filters = [this] (const auto& read) {
        return !passes_all_basic_filters(read);
    };
    
    if (!basic_filters_.empty()) {
        if (context_filters_.empty()) {
            last = std::partition(first, last, passes_basic_filters);
        } else {
            last = std::stable_partition(first, last, passes_basic_filters);
        }
    }
    
    std::for_each(cbegin(context_filters_), cend(context_filters_),
                  [first, &last] (const auto& filter) {
                      last = filter->partition(first, last);
                  });
    
    return last;
}

template <typename BidirIt>
BidirIt ReadFilter<BidirIt>::remove(BidirIt first, BidirIt last,
                                    FilterCountMap& filter_counts) const
{
    if (num_filters() == 0) return last;
    
    if (first == last) {
        filter_counts.reserve(num_filters());
        for (const auto& filter : basic_filters_) {
            filter_counts.emplace(filter->name(), 0);
        }
        for (const auto& filter : context_filters_) {
            filter_counts.emplace(filter->name(), 0);
        }
        return last;
    }
    
    if (!basic_filters_.empty()) {
        std::vector<std::size_t> flat_counts(basic_filters_.size(), 0);
        
        last = std::remove_if(first, last,
                              [this, &flat_counts] (const auto& read) {
                                  const auto it = find_failing_basic_filter(read);
                                  
                                  if (it != std::cend(basic_filters_)) {
                                      ++flat_counts[std::distance(std::cbegin(basic_filters_), it)];
                                      return true;
                                  }
                                  
                                  return false;
                              });
        
        filter_counts.reserve(num_filters());
        
        std::transform(std::cbegin(basic_filters_), std::cend(basic_filters_), std::cbegin(flat_counts),
                       std::inserter(filter_counts, std::begin(filter_counts)),
                       [] (const auto& filter, const auto count) {
                           return std::make_pair(filter->name(), count);
                       });
    }
    
    std::for_each(cbegin(context_filters_), cend(context_filters_),
                  [first, &last, &filter_counts] (const auto& filter) {
                      const auto it = filter->remove(first, last);
                      
                      filter_counts.emplace(filter->name(), std::distance(it, last));
                      
                      last = it;
                  });
    
    return last;
}

template <typename BidirIt>
BidirIt ReadFilter<BidirIt>::partition(BidirIt first, BidirIt last,
                                       FilterCountMap& filter_counts) const
{
    if (num_filters() == 0) return last;
    
    if (first == last) {
        filter_counts.reserve(num_filters());
        for (const auto& filter : basic_filters_) {
            filter_counts.emplace(filter->name(), 0);
        }
        for (const auto& filter : context_filters_) {
            filter_counts.emplace(filter->name(), 0);
        }
        return last;
    }
    
    if (!basic_filters_.empty()) {
        std::vector<std::size_t> flat_counts(basic_filters_.size(), 0);
        
        last = std::stable_partition(first, last,
                      [this, &flat_counts] (const auto& read) {
                          const auto it = find_failing_basic_filter(read);
                          
                          if (it != std::cend(basic_filters_)) {
                              ++flat_counts[std::distance(std::cbegin(basic_filters_), it)];
                              return true;
                          }
                          
                          return false;
                      });
        
        filter_counts.reserve(num_filters());
        
        std::transform(std::cbegin(basic_filters_), std::cend(basic_filters_), std::cbegin(flat_counts),
                       std::inserter(filter_counts, std::begin(filter_counts)),
                       [] (const auto& filter, const auto count) {
                           return std::make_pair(filter->name(), count);
                       });
    }
    
    std::for_each(cbegin(context_filters_), cend(context_filters_),
                  [first, &last, &filter_counts] (const auto& filter) {
                      const auto it = filter->partition(first, last);
                      
                      filter_counts.emplace(filter->name(), std::distance(it, last));
                      
                      last = it;
                  });
    
    return last;
}

// private member methods

template <typename BidirIt>
bool ReadFilter<BidirIt>::passes_all_basic_filters(const AlignedRead& read) const noexcept
{
    return std::all_of(std::cbegin(basic_filters_), std::cend(basic_filters_),
                       [&read] (const auto& filter) { return (*filter)(read); });
}

template <typename BidirIt>
auto ReadFilter<BidirIt>::find_failing_basic_filter(const AlignedRead& read) const noexcept
{
    return std::find_if_not(std::cbegin(basic_filters_), std::cend(basic_filters_),
                            [&read] (const auto& filter) { return (*filter)(read); });
}

// non-member methods

template <typename ReadFilter>
auto remove(MappableFlatMultiSet<AlignedRead>& reads, const ReadFilter& filter)
{
    return filter.remove(std::begin(reads), std::end(reads));
}

template <typename ReadFilter>
auto remove(MappableFlatMultiSet<AlignedRead>& reads, const ReadFilter& filter,
            typename ReadFilter::FilterCountMap& filter_counts)
{
    return filter.remove(std::begin(reads), std::end(reads), filter_counts);
}

template <typename KeyType>
using FilterPointMap = std::unordered_map<KeyType, typename MappableMap<KeyType, AlignedRead>::mapped_type::iterator>;

template <typename KeyType, typename ReadFilter>
FilterPointMap<KeyType>
filter(MappableMap<KeyType, AlignedRead>& reads, const ReadFilter& f)
{
    FilterPointMap<KeyType> result {reads.size()};
    
    if (true) {
        typename ReadFilter::FilterCountMap filter_counts {};
        
        for (auto& p : reads) {
            result.emplace(p.first, remove(p.second, f, filter_counts));
        }
        
        Logging::DebugLogger log {};
        
        for (const auto& p : filter_counts) {
            stream(log) << p.second << " reads were removed by the " << p.first << " filter";
        }
    } else {
        for (auto& p : reads) {
            result.emplace(p.first, remove(p.second, f));
        }
    }
    
    return result;
}

template <typename KeyType>
std::size_t erase_filtered_reads(MappableMap<KeyType, AlignedRead>& reads,
                                 const FilterPointMap<KeyType>& filter_points)
{
    std::size_t result {0};
    
    for (auto& p : reads) {
        const auto it = filter_points.at(p.first);
        result += std::distance(it, std::end(p.second));
        p.second.erase(it, std::end(p.second));
    }
    
    return result;
}

} // namespace Octopus

#endif /* defined(__Octopus__read_filter__) */
