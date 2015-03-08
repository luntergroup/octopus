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

using std::cbegin;
using std::cend;

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
    template <typename OutputIterator1, typename OutputIterator2>
    void filter_reads(BidirectionalIterator first, BidirectionalIterator last,
                      OutputIterator1 good_reads, OutputIterator2 bad_reads) const;
    
private:
    std::vector<ContextFreeFilter> context_free_filters_;
    std::vector<ContextBasedFilter> context_based_filters_;
    
    bool filter_read(const AlignedRead& the_read, BidirectionalIterator first,
                     BidirectionalIterator previous) const;
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

template <typename BidirectionalIterator>
template <typename OutputIterator1, typename OutputIterator2>
void ReadFilter<BidirectionalIterator>::filter_reads(BidirectionalIterator first,
                                                     BidirectionalIterator last,
                                                     OutputIterator1 good_reads,
                                                     OutputIterator2 bad_reads) const
{
    BidirectionalIterator previous {first};
    std::partition_copy(first, last, good_reads, bad_reads,
                        [this, first, &previous] (const AlignedRead& the_read) {
                            return filter_read(the_read, first, (previous != first) ? previous++ :
                                               previous);
                        });
}

template <typename BidirectionalIterator>
bool ReadFilter<BidirectionalIterator>::filter_read(const AlignedRead& the_read,
                                                    BidirectionalIterator first,
                                                    BidirectionalIterator previous) const
{
    return std::all_of(cbegin(context_free_filters_), cend(context_free_filters_),
                       [&the_read] (const auto& filter) {
                           return filter(the_read);
    }) && std::all_of(cbegin(context_based_filters_), cend(context_based_filters_),
                      [&the_read, first, previous] (const auto& filter) {
                          return filter(the_read, first, previous);
    });
}

#endif /* defined(__Octopus__read_filter__) */
