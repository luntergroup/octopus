//
//  read_filter.h
//  Octopus
//
//  Created by Daniel Cooke on 06/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__read_filter__
#define __Octopus__read_filter__

#include <functional>
#include <vector>
#include <algorithm> // std::all_of, std::partition_copy
#include <iterator>  // std::cbegin etc
#include <type_traits>

#include "context_back_insert_iterator.h"

template <typename BidirectionalIterator>
class ReadFilter
{
public:
    using ContextFreeFilter  = std::function<bool(const AlignedRead&)>;
    using ContextBasedFilter = std::function<bool(const AlignedRead&, BidirectionalIterator,
                                                  BidirectionalIterator)>;
    
    ReadFilter()  = default;
    ~ReadFilter() = default;
    
    ReadFilter(const ReadFilter&)            = default;
    ReadFilter& operator=(const ReadFilter&) = default;
    ReadFilter(ReadFilter&&)                 = default;
    ReadFilter& operator=(ReadFilter&&)      = default;
    
    void register_filter(ContextFreeFilter a_filter);
    void register_filter(ContextBasedFilter a_filter);
    unsigned num_filters() const noexcept;
    template <typename InputIterator, typename OutputIterator1, typename OutputIterator2>
    std::pair<OutputIterator1, OutputIterator2>
    filter_reads(InputIterator first, InputIterator last,
                 OutputIterator1 good_reads, OutputIterator2 bad_reads) const;
    
private:
    std::vector<ContextFreeFilter> context_free_filters_;
    std::vector<ContextBasedFilter> context_based_filters_;
    
    bool filter_read(const AlignedRead& the_read, BidirectionalIterator first_good,
                     BidirectionalIterator previous_good) const;
};

template <typename BidirectionalIterator>
void ReadFilter<BidirectionalIterator>::register_filter(ContextFreeFilter a_filter)
{
    context_free_filters_.emplace_back(std::move(a_filter));
}

template <typename BidirectionalIterator>
void ReadFilter<BidirectionalIterator>::register_filter(ContextBasedFilter a_filter)
{
    context_based_filters_.emplace_back(std::move(a_filter));
}

template <typename BidirectionalIterator>
unsigned ReadFilter<BidirectionalIterator>::num_filters() const noexcept
{
    return static_cast<unsigned>(context_free_filters_.size() + context_based_filters_.size());
}

template <typename T>
inline
typename std::enable_if<std::is_void<typename T::value_type>::value,
    typename T::container_type::const_iterator>::type
get_first(T first, T last)
{
    return last.begin();
}

template <typename T>
inline
typename std::enable_if<!std::is_void<typename T::value_type>::value, T>::type
get_first(T first, T last)
{
    return first;
}

template <typename T>
inline
typename std::enable_if<std::is_void<typename T::value_type>::value,
    typename T::container_type::const_iterator>::type
get_last(T first, T last)
{
    return (last.begin() != last.end()) ? std::prev(last.end()) : last.begin();
}

template <typename T>
inline
typename std::enable_if<!std::is_void<typename T::value_type>::value, T>::type
get_last(T first, T last)
{
    return (first != last) ? std::prev(last) : last;
}

template <typename BidirectionalIterator>
template <typename InputIterator, typename OutputIterator1, typename OutputIterator2>
std::pair<OutputIterator1, OutputIterator2>
ReadFilter<BidirectionalIterator>::filter_reads(InputIterator first, InputIterator last,
                                                 OutputIterator1 good_reads,
                                                 OutputIterator2 bad_reads) const
{
    auto good_reads_last = good_reads;
    return std::partition_copy(first, last, good_reads, bad_reads,
        [this, good_reads, &good_reads_last] (const AlignedRead& the_read) {
        if (filter_read(the_read, get_first(good_reads, good_reads_last),
                        get_last(good_reads, good_reads_last))) {
            ++good_reads_last;
            return true;
        }
        return false;
    });
}

template <typename BidirectionalIterator>
bool ReadFilter<BidirectionalIterator>::filter_read(const AlignedRead& the_read,
                                                    BidirectionalIterator first_good,
                                                    BidirectionalIterator previous_good) const
{
    return std::all_of(std::cbegin(context_free_filters_), std::cend(context_free_filters_),
                       [&the_read] (const auto& filter) {
                           return filter(the_read);
    }) && std::all_of(std::cbegin(context_based_filters_), std::cend(context_based_filters_),
                      [&the_read, first_good, previous_good] (const auto& filter) {
                          return filter(the_read, first_good, previous_good);
    });
}

#endif /* defined(__Octopus__read_filter__) */
