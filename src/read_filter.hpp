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
#include <functional>  // std::function
#include <algorithm>   // std::all_of, std::partition_copy
#include <iterator>    // std::prev, std::cbegin, std::cend
#include <type_traits> // std::enable_if, std::is_void
#include <utility>     // std::move

namespace Octopus {

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
    
    void register_filter(ContextFreeFilter filter);
    void register_filter(ContextBasedFilter filter);
    
    unsigned num_filters() const noexcept;
    
    template <typename InputIterator, typename OutputIterator1, typename OutputIterator2>
    std::pair<OutputIterator1, OutputIterator2>
    filter_reads(InputIterator first, InputIterator last,
                 OutputIterator1 good_reads, OutputIterator2 bad_reads) const;
    
private:
    std::vector<ContextFreeFilter> noncontext_filters_;
    std::vector<ContextBasedFilter> context_filters_;
    
    bool filter_read(const AlignedRead& read, BidirectionalIterator first_good,
                     BidirectionalIterator previous_good) const;
};

template <typename BidirectionalIterator>
void ReadFilter<BidirectionalIterator>::register_filter(ContextFreeFilter filter)
{
    noncontext_filters_.emplace_back(std::move(filter));
}

template <typename BidirectionalIterator>
void ReadFilter<BidirectionalIterator>::register_filter(ContextBasedFilter filter)
{
    context_filters_.emplace_back(std::move(filter));
}

template <typename BidirectionalIterator>
unsigned ReadFilter<BidirectionalIterator>::num_filters() const noexcept
{
    return static_cast<unsigned>(noncontext_filters_.size() + context_filters_.size());
}

namespace detail
{
    template <typename T>
    inline
    typename std::enable_if<std::is_void<typename T::value_type>::value, typename T::container_type::const_iterator>::type
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
    typename std::enable_if<std::is_void<typename T::value_type>::value, typename T::container_type::const_iterator>::type
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
} // namespace detail

template <typename BidirectionalIterator>
template <typename InputIterator, typename OutputIterator1, typename OutputIterator2>
std::pair<OutputIterator1, OutputIterator2>
ReadFilter<BidirectionalIterator>::filter_reads(InputIterator first, InputIterator last,
                                                OutputIterator1 good_reads, OutputIterator2 bad_reads) const
{
    auto good_reads_last = good_reads;
    return std::partition_copy(first, last, good_reads, bad_reads,
        [this, good_reads, &good_reads_last] (const AlignedRead& read) {
            if (filter_read(read, detail::get_first(good_reads, good_reads_last),
                            detail::get_last(good_reads, good_reads_last))) {
            ++good_reads_last;
            return true;
        }
        return false;
    });
}

template <typename Iterator>
bool
ReadFilter<Iterator>::filter_read(const AlignedRead& read, Iterator first_good, Iterator prev_good) const
{
    using std::cbegin; using std::cend; using std::all_of;
    
    return all_of(cbegin(noncontext_filters_), cend(noncontext_filters_),
                       [&read] (const auto& filter) { return filter(read); })
        && all_of(cbegin(context_filters_), cend(context_filters_),
                  [&read, first_good, prev_good] (const auto& filter) { return filter(read, first_good, prev_good); });
}

} // namespace Octopus

#endif /* defined(__Octopus__read_filter__) */
