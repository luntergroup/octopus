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
#include <functional>
#include <algorithm>
#include <iterator>
#include <type_traits>
#include <utility>
#include <unordered_map>

#include <boost/optional.hpp>

#include "mappable_flat_multi_set.hpp"
#include "mappable_map.hpp"
#include "logging.hpp"

namespace Octopus
{
/*
 ReadFilter stores a collection of filter functions, which can either be non-context-based (can be 
 applied to a single read), or context-based (conditional on what has already been filtered).
 
 This is a template class as the type of iterator used for context-based filteration needs to be 
 known at compile time. Essentially the class needs to know what container it is going to be operating on.
 */
template <typename BidirIt>
class ReadFilter
{
public:
    using Iterator = BidirIt;
    
    using ContextFreeFilter  = std::function<bool(const AlignedRead&)>;
    using ContextBasedFilter = std::function<Iterator(Iterator, Iterator)>;
    
    ReadFilter()  = default;
    ReadFilter(Logging::InfoLogger log) : log_ {log} {};
    ~ReadFilter() = default;
    
    ReadFilter(const ReadFilter&)            = default;
    ReadFilter& operator=(const ReadFilter&) = default;
    ReadFilter(ReadFilter&&)                 = default;
    ReadFilter& operator=(ReadFilter&&)      = default;
    
    void register_filter(ContextFreeFilter filter);
    void register_filter(ContextBasedFilter filter);
    
    unsigned num_filters() const noexcept;
    
    Iterator filter(Iterator first, Iterator last) const;
    
private:
    std::vector<ContextFreeFilter> noncontext_filters_;
    std::vector<ContextBasedFilter> context_filters_;
    
    boost::optional<Logging::InfoLogger> log_;
    
    bool passes_all_noncontext_filters(const AlignedRead& read) const;
    
    Iterator remove(Iterator first, Iterator last) const;
    Iterator partition(Iterator first, Iterator last) const;
};

template <typename BidirIt>
void ReadFilter<BidirIt>::register_filter(ContextFreeFilter filter)
{
    noncontext_filters_.emplace_back(std::move(filter));
}

template <typename BidirIt>
void ReadFilter<BidirIt>::register_filter(ContextBasedFilter filter)
{
    context_filters_.emplace_back(std::move(filter));
}

template <typename BidirIt>
unsigned ReadFilter<BidirIt>::num_filters() const noexcept
{
    return static_cast<unsigned>(noncontext_filters_.size() + context_filters_.size());
}

template <typename BidirIt>
BidirIt ReadFilter<BidirIt>::filter(BidirIt first, BidirIt last) const
{
    if (log_) {
        return partition(first, last);
    }
    
    return remove(first, last);
}

// private member methods

template <typename Iterator>
bool ReadFilter<Iterator>::passes_all_noncontext_filters(const AlignedRead& read) const
{
    return std::all_of(std::cbegin(noncontext_filters_), std::cend(noncontext_filters_),
                       [&read] (const auto& filter) { return filter(read); });
}

template <typename BidirIt>
BidirIt ReadFilter<BidirIt>::remove(BidirIt first, BidirIt last) const
{
    using std::all_of; using std::cbegin; using std::cend;
    
    auto it = std::remove_if(first, last,
                             [this] (const auto& read) {
                                 return !passes_all_noncontext_filters(read);
                             });
    
    std::for_each(cbegin(context_filters_), cend(context_filters_),
                  [first, &it] (const auto& filter) { it = filter(first, it); });
    
    return it;
}

template <typename BidirIt>
BidirIt ReadFilter<BidirIt>::partition(BidirIt first, BidirIt last) const
{
    using std::all_of; using std::cbegin; using std::cend;
    
    auto it = std::stable_partition(first, last,
                                    [this] (const auto& read) {
                                        return !passes_all_noncontext_filters(read);
                                    });
    
    std::unordered_map<std::string, std::size_t> failed_filter_counts {noncontext_filters_.size()};
    
//    for (const auto& filter : noncontext_filters_) {
//        failed_filter_counts.emplace(filter.name(), 0);
//    }
//    
//    std::for_each(it, last,
//                  [this] (const auto& read) {
//                      for (const auto& filter : noncontext_filters_) {
//                          if (filter(read)) ++failed_filter_counts.at(filter.name());
//                      }
//                  });
//    
//    for (const auto& p : failed_filter_counts) {
//        if (p.second > 0) {
//            stream(*log_) << p.second << " reads failed the " << p.first << " filter";
//        }
//    }
    
    std::for_each(cbegin(context_filters_), cend(context_filters_),
                  [first, &it] (const auto& filter) { it = filter(first, it); });
    
    return it;
}

// non-member methods

template <typename ReadFilter>
auto filter(MappableFlatMultiSet<AlignedRead>& reads, const ReadFilter& filter)
{
    return filter.filter(std::begin(reads), std::end(reads));
}

template <typename KeyType>
using PartitionPointMap = std::unordered_map<KeyType, typename MappableMap<KeyType, AlignedRead>::mapped_type::iterator>;

template <typename KeyType, typename ReadFilter>
PartitionPointMap<KeyType>
filter(MappableMap<KeyType, AlignedRead>& reads, const ReadFilter& f)
{
    PartitionPointMap<KeyType> result {reads.size()};
    
    for (auto& p : reads) {
        result.emplace(p.first, filter(p.second, f));
    }
    
    return result;
}

template <typename KeyType>
std::size_t erase_filtered_reads(MappableMap<KeyType, AlignedRead>& reads,
                                 const PartitionPointMap<KeyType>& partition_points)
{
    std::size_t result {0};
    
    for (auto& p : reads) {
        const auto it = partition_points.at(p.first);
        result += std::distance(it, std::end(p.second));
        p.second.erase(it, std::end(p.second));
    }
    
    return result;
}

} // namespace Octopus

#endif /* defined(__Octopus__read_filter__) */
