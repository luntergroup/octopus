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

#include <iostream> // DEBUG

namespace Octopus {

/*
 ReadFilter stores a collection of filter functions, which can either be non-context-based (can be 
 applied to a single read), or context-based (conditional on what has already been filtered).
 
 This is a template class as the type of iterator used for context-based filteration needs to be 
 known at compile time. Essentially the class needs to know what container it is going to be operating on.
 */
template <typename BidirectionalIterator>
class ReadFilter
{
public:
    using Iterator = BidirectionalIterator;
    
    using ContextFreeFilter  = std::function<bool(const AlignedRead&)>;
    using ContextBasedFilter = std::function<bool(const AlignedRead&, Iterator, Iterator)>;
    
    ReadFilter()  = default;
    ~ReadFilter() = default;
    
    ReadFilter(const ReadFilter&)            = default;
    ReadFilter& operator=(const ReadFilter&) = default;
    ReadFilter(ReadFilter&&)                 = default;
    ReadFilter& operator=(ReadFilter&&)      = default;
    
    void register_filter(ContextFreeFilter filter);
    void register_filter(ContextBasedFilter filter);
    
    unsigned num_filters() const noexcept;
    
    // removes all failing reads, which may not be preserved
    Iterator filter(Iterator first, Iterator last) const;
    
    // partitions in-place, may use aditional memory if available to improve performance
    //Iterator partition(Iterator first, Iterator last) const;
    
    // copy partitions
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
    using IsContextIterator = std::is_void<typename T::value_type>;
    
    template <typename T>
    auto get_first(T first, T last, std::true_type)
    {
        return std::begin(last);
    }
    
    template <typename T>
    auto get_first(T first, T last, std::false_type)
    {
        return first;
    }
    
    template <typename T>
    auto get_first(T first, T last)
    {
        return get_first(first, last, IsContextIterator<T> {});
    }
    
    template <typename T>
    auto get_last(T first, T last, std::true_type)
    {
        using std::begin; using std::end;
        return (begin(last) != end(last)) ? std::prev(end(last)) : begin(last);
    }
    
    template <typename T>
    auto get_last(T first, T last, std::false_type)
    {
        return (first != last) ? std::prev(last) : last;
    }
    
    template <typename T>
    auto get_last(T first, T last)
    {
        return get_last(first, last, IsContextIterator<T> {});
    }
} // namespace detail

template <typename BidirectionalIterator>
BidirectionalIterator
ReadFilter<BidirectionalIterator>::filter(BidirectionalIterator first, BidirectionalIterator last) const
{
    using std::remove_if; using std::all_of; using std::cbegin; using std::cend;
    
    auto it = remove_if(first, last, [this] (const auto& read) {
        return all_of(cbegin(noncontext_filters_), cend(noncontext_filters_),
                      [&read] (const auto& filter) { return filter(read); });
    });
    
    auto prev = first;
    
    return remove_if(first, it, [this, first, &prev] (const auto& read) {
        auto all_good = all_of(cbegin(context_filters_), cend(context_filters_),
                               [&read, first, prev] (const auto& filter) {
                                   return filter(read, first, prev);
                               });
        
        if (all_good) {
            ++prev;
        }
        
        return all_good;
    });
}

//template <typename BidirectionalIterator>
//BidirectionalIterator
//ReadFilter<BidirectionalIterator>::partition(BidirectionalIterator first, BidirectionalIterator last) const
//{
//    using std::stable_partition; using std::all_of; using std::cbegin; using std::cend;
//    
//    auto it = stable_partition(first, last, [this] (const auto& read) {
//        return all_of(cbegin(noncontext_filters_), cend(noncontext_filters_),
//                      [&read] (const auto& filter) { return filter(read); });
//    });
//    
//    // now we have potentially it - begin passes
//    
//    std::vector<typename BidirectionalIterator::value_type> buffer(std::distance(first, it));
//    std::move(first, it, begin(buffer));
//    
//    auto prev = first;
//    
//    std::partition_copy(std::make_move_iterator(begin(buffer)),
//                        std::make_move_iterator(end(buffer)),
//                        first, std::reverse_iterator<BidirectionalIterator> {it},
//                        [this, first, &prev] (const auto& read) {
//                            auto all_good = all_of(cbegin(context_filters_), cend(context_filters_),
//                                                   [&read, first, prev] (const auto& filter) {
//                                                       return filter(read, detail::get_first(first, prev),
//                                                                     detail::get_last(first, prev));
//                                                   });
//                            
//                            if (all_good) ++prev;
//                            
//                            return all_good;
//                        });
//    
//    return prev;
//}

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
bool ReadFilter<Iterator>::filter_read(const AlignedRead& read, Iterator first_good, Iterator prev_good) const
{
    using std::cbegin; using std::cend; using std::all_of;
    
    return all_of(cbegin(noncontext_filters_), cend(noncontext_filters_),
                       [&read] (const auto& filter) { return filter(read); })
        && all_of(cbegin(context_filters_), cend(context_filters_),
                  [&read, first_good, prev_good] (const auto& filter) {
                      return filter(read, first_good, prev_good);
                  });
}

} // namespace Octopus

#endif /* defined(__Octopus__read_filter__) */
