//
//  read_filterer.hpp
//  Octopus
//
//  Created by Daniel Cooke on 06/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__read_filterer__
#define __Octopus__read_filterer__

#include <vector>
#include <unordered_map>
#include <memory>
#include <algorithm>
#include <iterator>
#include <type_traits>
#include <utility>

#include <boost/optional.hpp>

#include "read_filter.hpp"

namespace octopus { namespace readpipe
{
/*
 ReadFilterer stores a collection of filter functions, which can either be non-context-based (
 applied to a single read), or context-based (applied to a range of reads).
 
 This is a template class as the type of iterator used for context-based filteration needs to be 
 known at compile time. The class needs to know what container it is going to be operating on.
 */
template <typename BidirIt>
class ReadFilterer
{
public:
    using ReadIterator = BidirIt;
    
    using BasicFilterPtr   = std::unique_ptr<BasicReadFilter>;
    using ContextFilterPtr = std::unique_ptr<ContextReadFilter<ReadIterator>>;
    
    using FilterCountMap = std::unordered_map<std::string, std::size_t>;
    
    ReadFilterer() = default;
    
    ReadFilterer(const ReadFilterer&)            = default;
    ReadFilterer& operator=(const ReadFilterer&) = default;
    ReadFilterer(ReadFilterer&&)                 = default;
    ReadFilterer& operator=(ReadFilterer&&)      = default;
    
    ~ReadFilterer() = default;
    
    void add(BasicFilterPtr filter);
    void add(ContextFilterPtr filter);
    
    unsigned num_filters() const noexcept;
    
    void shrink_to_fit() noexcept; // Just removes extra capcity for filters
    
    // Like std::remove
    BidirIt remove(ReadIterator first, ReadIterator last) const;
    BidirIt remove(ReadIterator first, ReadIterator last, FilterCountMap& filter_counts) const;
    
    // Like std::stable_partition
    BidirIt partition(ReadIterator first, ReadIterator last) const;
    BidirIt partition(ReadIterator first, ReadIterator last, FilterCountMap& filter_counts) const;
    
private:
    std::vector<BasicFilterPtr> basic_filters_;
    std::vector<ContextFilterPtr> context_filters_;
    
    bool passes_all_basic_filters(const AlignedRead& read) const noexcept;
    auto find_failing_basic_filter(const AlignedRead& read) const noexcept;
};

template <typename BidirIt>
void ReadFilterer<BidirIt>::add(BasicFilterPtr filter)
{
    basic_filters_.emplace_back(std::move(filter));
}

template <typename BidirIt>
void ReadFilterer<BidirIt>::add(ContextFilterPtr filter)
{
    context_filters_.emplace_back(std::move(filter));
}

template <typename BidirIt>
unsigned ReadFilterer<BidirIt>::num_filters() const noexcept
{
    return static_cast<unsigned>(basic_filters_.size() + context_filters_.size());
}

template <typename BidirIt>
void ReadFilterer<BidirIt>::shrink_to_fit() noexcept
{
    basic_filters_.shrink_to_fit();
    context_filters_.shrink_to_fit();
}

template <typename BidirIt>
BidirIt ReadFilterer<BidirIt>::remove(BidirIt first, BidirIt last) const
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
BidirIt ReadFilterer<BidirIt>::partition(BidirIt first, BidirIt last) const
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
BidirIt ReadFilterer<BidirIt>::remove(BidirIt first, BidirIt last,
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
BidirIt ReadFilterer<BidirIt>::partition(BidirIt first, BidirIt last,
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
bool ReadFilterer<BidirIt>::passes_all_basic_filters(const AlignedRead& read) const noexcept
{
    return std::all_of(std::cbegin(basic_filters_), std::cend(basic_filters_),
                       [&read] (const auto& filter) { return (*filter)(read); });
}

template <typename BidirIt>
auto ReadFilterer<BidirIt>::find_failing_basic_filter(const AlignedRead& read) const noexcept
{
    return std::find_if_not(std::cbegin(basic_filters_), std::cend(basic_filters_),
                            [&read] (const auto& filter) { return (*filter)(read); });
}

// non-member methods

template <typename Container>
using ReadFiltererTp = ReadFilterer<typename Container::iterator>;

template <typename Container, typename ReadFilterer>
auto remove(Container& reads, const ReadFilterer& filter)
{
    return filter.remove(std::begin(reads), std::end(reads));
}

template <typename Container, typename ReadFilterer>
auto remove(Container& reads, const ReadFilterer& filter,
            typename ReadFilterer::FilterCountMap& filter_counts)
{
    return filter.remove(std::begin(reads), std::end(reads), filter_counts);
}

template <typename Map>
using FilterPointMap = std::unordered_map<typename Map::key_type, typename Map::mapped_type::iterator>;

template <typename ReadFilterer>
using FilterCountMap = typename ReadFilterer::FilterCountMap;

template <typename S, typename ReadFilterer>
using SampleFilterCountMap = std::unordered_map<S, FilterCountMap<ReadFilterer>>;

template <typename Map, typename ReadFilterer>
using OptionalFilterCountMap = boost::optional<SampleFilterCountMap<typename Map::key_type, ReadFilterer>&>;

template <typename Map, typename ReadFilterer>
FilterPointMap<Map>
filter(Map& reads, const ReadFilterer& f, OptionalFilterCountMap<Map, ReadFilterer> filter_counts = boost::none)
{
    FilterPointMap<Map> result {reads.size()};
    
    for (auto& p : reads) {
        if (filter_counts && filter_counts->count(p.first) == 1) {
            result.emplace(p.first, remove(p.second, f, filter_counts->at(p.first)));
        } else {
            result.emplace(p.first, remove(p.second, f));
        }
    }
    
    return result;
}

template <typename Map>
std::size_t erase_filtered_reads(Map& reads, const FilterPointMap<Map>& filter_points)
{
    std::size_t result {0};
    
    for (auto& p : reads) {
        const auto it = filter_points.at(p.first);
        result += std::distance(it, std::end(p.second));
        p.second.erase(it, std::end(p.second));
    }
    
    return result;
}
} // namespace readpipe
} // namespace octopus

#endif /* defined(__Octopus__read_filterer__) */
