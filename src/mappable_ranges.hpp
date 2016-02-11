//
//  mappable_ranges.h
//  Octopus
//
//  Created by Daniel Cooke on 20/07/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_mappable_ranges_h
#define Octopus_mappable_ranges_h

#include <iterator>
#include <cstddef>
#include <type_traits>

#include <boost/iterator/filter_iterator.hpp>
#include <boost/range/iterator_range_core.hpp>

/*
 An ordered collection of MappableType elements, X, is:
 - ForwardSorted         iff i <= j -> get_region(X[i]) <= get_region(X[j])
 - BidirectionallySorted iff X is ForwardSorted AND i <= j -> get_end(X[i]) <= get_end(X[j])
 - Unsorted              iff X is not ForwardSorted
 */
enum class MappableRangeOrder { ForwardSorted, BidirectionallySorted, Unsorted };

namespace detail
{
    template <typename MappableType>
    class IsOverlapped
    {
    public:
        IsOverlapped() = delete;
        template <typename MappableType2>
        IsOverlapped(const MappableType2& mappable) : region_ {mapped_region(mappable)} {}
        
        bool operator()(const MappableType& mappable) { return overlaps(mappable, region_); }
        
    private:
        RegionType<MappableType> region_;
    };
} // namespace detail

template <typename Iterator>
using OverlapIterator = boost::filter_iterator<detail::IsOverlapped<typename std::iterator_traits<Iterator>::value_type>, Iterator>;

template <typename Iterator>
inline bool operator==(Iterator lhs, OverlapIterator<Iterator> rhs) noexcept
{
    return lhs == rhs.base();
}
template <typename Iterator>
inline bool operator!=(Iterator lhs, OverlapIterator<Iterator> rhs) noexcept
{
    return !operator==(lhs, rhs);
}
template <typename Iterator>
inline bool operator==(OverlapIterator<Iterator> lhs, Iterator rhs) noexcept
{
    return operator==(rhs, lhs);
}
template <typename Iterator>
inline bool operator!=(OverlapIterator<Iterator> lhs, Iterator rhs) noexcept
{
    return !operator==(lhs, rhs);
}

template <typename Iterator>
using OverlapRange = boost::iterator_range<OverlapIterator<Iterator>>;

template <typename Iterator>
inline
boost::iterator_range<Iterator> bases(const OverlapRange<Iterator>& overlap_range)
{
    return boost::make_iterator_range(overlap_range.begin().base(), overlap_range.end().base());
}

template <typename Iterator>
inline
size_t size(const OverlapRange<Iterator>& range, MappableRangeOrder order = MappableRangeOrder::ForwardSorted)
{
    if (order == MappableRangeOrder::BidirectionallySorted) {
        return std::distance(range.begin().base(), range.end().base());
    } else {
        return std::distance(range.begin(), range.end());
    }
}

template <typename Iterator>
inline
bool empty(const OverlapRange<Iterator>& range)
{
    return range.empty();
}

template <typename Iterator, typename MappableType>
inline
OverlapRange<Iterator>
make_overlap_range(Iterator first, Iterator last, const MappableType& mappable)
{
    using boost::make_iterator_range; using boost::make_filter_iterator; using detail::IsOverlapped;
    
    using MappableType2 = typename std::iterator_traits<Iterator>::value_type;
    
    return make_iterator_range(
        make_filter_iterator<IsOverlapped<MappableType2>>(IsOverlapped<MappableType2>(mappable), first, last),
        make_filter_iterator<IsOverlapped<MappableType2>>(IsOverlapped<MappableType2>(mappable), last, last)
    );
}

namespace detail
{
    template <typename MappableType>
    class IsContained
    {
    public:
        IsContained() = delete;
        template <typename MappableType_>
        IsContained(const MappableType_& mappable) : region_ {mapped_region(mappable)} {}
        
        bool operator()(const MappableType& mappable) { return contains(region_, mappable); }
        
    private:
        RegionType<MappableType> region_;
    };
} // namespace detail

template <typename Iterator>
using ContainedIterator = boost::filter_iterator<detail::IsContained<typename std::iterator_traits<Iterator>::value_type>, Iterator>;

template <typename Iterator>
inline bool operator==(Iterator lhs, ContainedIterator<Iterator> rhs) noexcept
{
    return lhs == rhs.base();
}
template <typename Iterator>
inline bool operator!=(Iterator lhs, ContainedIterator<Iterator> rhs) noexcept
{
    return !operator==(lhs, rhs);
}
template <typename Iterator>
inline bool operator==(ContainedIterator<Iterator> lhs, Iterator rhs) noexcept
{
    return operator==(rhs, lhs);
}
template <typename Iterator>
inline bool operator!=(ContainedIterator<Iterator> lhs, Iterator rhs) noexcept
{
    return !operator==(lhs, rhs);
}

template <typename Iterator>
using ContainedRange = boost::iterator_range<ContainedIterator<Iterator>>;

template <typename Iterator>
inline
boost::iterator_range<Iterator> bases(const ContainedRange<Iterator>& contained_range)
{
    return boost::make_iterator_range(contained_range.begin().base(), contained_range.end().base());
}

template <typename Iterator>
inline
size_t size(const ContainedRange<Iterator>& range, MappableRangeOrder order = MappableRangeOrder::ForwardSorted)
{
    if (order == MappableRangeOrder::BidirectionallySorted) {
        return std::distance(range.begin().base(), range.end().base());
    } else {
        return std::distance(range.begin(), range.end());
    }
}

template <typename Iterator>
inline
bool empty(const ContainedRange<Iterator>& range)
{
    return range.empty();
}

template <typename Iterator, typename MappableType>
inline
ContainedRange<Iterator>
make_contained_range(Iterator first, Iterator last, const MappableType& mappable)
{
    using boost::make_iterator_range; using boost::make_filter_iterator; using detail::IsContained;
    
    using MappableType2 = typename std::iterator_traits<Iterator>::value_type;
    
    return make_iterator_range(
        make_filter_iterator<IsContained<MappableType2>>(IsContained<MappableType2>(mappable), first, last),
        make_filter_iterator<IsContained<MappableType2>>(IsContained<MappableType2>(mappable), last, last)
    );
}

namespace detail
{
    template <typename MappableType>
    class IsShared
    {
    public:
        IsShared() = delete;
        template <typename MappableType1_, typename MappableType2_>
        IsShared(const MappableType1_& lhs, MappableType2_ rhs)
        : lhs_ {get_region(lhs)}, rhs_ {get_region(rhs)} {}
        
        bool operator()(const MappableType& mappable)
        {
            return overlaps(lhs_, mappable) && overlaps(mappable, rhs_);
        }
        
    private:
        RegionType<MappableType> lhs_, rhs_;
    };
} // namespace detail

template <typename Iterator>
using SharedIterator = boost::filter_iterator<detail::IsShared<typename std::iterator_traits<Iterator>::value_type>, Iterator>;

template <typename Iterator>
inline bool operator==(Iterator lhs, SharedIterator<Iterator> rhs) noexcept
{
    return lhs == rhs.base();
}
template <typename Iterator>
inline bool operator!=(Iterator lhs, SharedIterator<Iterator> rhs) noexcept
{
    return !operator==(lhs, rhs);
}
template <typename Iterator>
inline bool operator==(SharedIterator<Iterator> lhs, Iterator rhs) noexcept
{
    return operator==(rhs, lhs);
}
template <typename Iterator>
inline bool operator!=(SharedIterator<Iterator> lhs, Iterator rhs) noexcept
{
    return !operator==(lhs, rhs);
}

template <typename Iterator>
using SharedRange = boost::iterator_range<SharedIterator<Iterator>>;

template <typename Iterator>
inline
boost::iterator_range<Iterator> bases(const SharedRange<Iterator>& shared_range)
{
    return boost::make_iterator_range(shared_range.begin().base(), shared_range.end().base());
}

template <typename Iterator>
inline
size_t size(const SharedRange<Iterator>& range, MappableRangeOrder order = MappableRangeOrder::ForwardSorted)
{
    if (order == MappableRangeOrder::BidirectionallySorted) {
        return std::distance(range.begin().base(), range.end().base());
    } else {
        return std::distance(range.begin(), range.end());
    }
}

template <typename Iterator>
inline
bool empty(const SharedRange<Iterator>& range)
{
    return range.empty();
}

template <typename Iterator, typename MappableType1, typename MappableType2>
inline
SharedRange<Iterator>
make_shared_range(Iterator first, Iterator last, const MappableType1& mappable1, const MappableType2& mappable2)
{
    using boost::make_iterator_range; using boost::make_filter_iterator; using detail::IsShared;
    
    using MappableType_ = typename std::iterator_traits<Iterator>::value_type;
    
    return make_iterator_range(
        make_filter_iterator<IsShared<MappableType_>>(IsShared<MappableType_>(mappable1, mappable2), first, last),
        make_filter_iterator<IsShared<MappableType_>>(IsShared<MappableType_>(mappable1, mappable2), last, last)
    );
}

#endif
