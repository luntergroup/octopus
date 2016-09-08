// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef mappable_ranges_hpp
#define mappable_ranges_hpp

#include <iterator>
#include <type_traits>
#include <cstddef>

#include <boost/iterator/filter_iterator.hpp>
#include <boost/range/iterator_range_core.hpp>

#include "mappable.hpp"

namespace octopus {

/*
 An ordered collection of MappableType elements, X, is:
 - ForwardSorted         iff i <= j -> mapped_region(X[i]) <= mapped_region(X[j])
 - BidirectionallySorted iff X is ForwardSorted AND i <= j -> end(X[i]) <= end(X[j])
 - Unsorted              iff X is not ForwardSorted
 */

struct ForwardSortedTag {};
struct BidirectionallySortedTag {};

namespace detail {

template <typename MappableType>
class IsOverlapped
{
public:
    IsOverlapped() = delete;
    
    template <typename MappableType2>
    IsOverlapped(const MappableType2& mappable) : region_ {mapped_region(mappable)} {}
    
    bool operator()(const MappableType& mappable) const
    {
        return overlaps(mappable, region_);
    }
    
private:
    RegionType<MappableType> region_;
};

} // namespace detail

template <typename Iterator>
using OverlapIterator = boost::filter_iterator<
    detail::IsOverlapped<typename std::iterator_traits<Iterator>::value_type>,
    Iterator
>;

template <typename Iterator>
using OverlapRange = boost::iterator_range<OverlapIterator<Iterator>>;

template <typename Iterator>
boost::iterator_range<Iterator> bases(const OverlapRange<Iterator>& range)
{
    return boost::make_iterator_range(range.begin().base(), range.end().base());
}

template <typename Iterator>
auto base_size(const OverlapRange<Iterator>& range)
{
    return static_cast<std::size_t>(std::distance(range.begin().base(), range.end().base()));
}

template <typename Iterator>
auto size(const OverlapRange<Iterator>& range, ForwardSortedTag)
{
    return static_cast<std::size_t>(std::distance(range.begin(), range.end()));
}

template <typename Iterator>
auto size(const OverlapRange<Iterator>& range, BidirectionallySortedTag)
{
    return base_size(range);
}

template <typename Iterator>
auto size(const OverlapRange<Iterator>& range)
{
    return size(range, ForwardSortedTag {});
}

template <typename Iterator>
bool empty(const OverlapRange<Iterator>& range)
{
    return range.empty();
}

template <typename Iterator, typename MappableType>
OverlapRange<Iterator> make_overlap_range(Iterator first, Iterator last, const MappableType& mappable)
{
    using boost::make_iterator_range; using boost::make_filter_iterator; using detail::IsOverlapped;
    using MappableType2 = typename std::iterator_traits<Iterator>::value_type;
    return make_iterator_range(
        make_filter_iterator<IsOverlapped<MappableType2>>(IsOverlapped<MappableType2>(mappable), first, last),
        make_filter_iterator<IsOverlapped<MappableType2>>(IsOverlapped<MappableType2>(mappable), last, last)
    );
}

namespace detail {

template <typename MappableType>
class IsContained
{
public:
    IsContained() = delete;
    
    template <typename MappableType_>
    IsContained(const MappableType_& mappable) : region_ {mapped_region(mappable)} {}
    
    bool operator()(const MappableType& mappable) const
    {
        return contains(region_, mappable);
    }
    
private:
    RegionType<MappableType> region_;
};

} // namespace detail

template <typename Iterator>
using ContainedIterator = boost::filter_iterator<detail::IsContained<
    typename std::iterator_traits<Iterator>::value_type>,
    Iterator
>;

template <typename Iterator>
using ContainedRange = boost::iterator_range<ContainedIterator<Iterator>>;

template <typename Iterator>
boost::iterator_range<Iterator> bases(const ContainedRange<Iterator>& range)
{
    return boost::make_iterator_range(range.begin().base(), range.end().base());
}

template <typename Iterator>
auto base_size(const ContainedRange<Iterator>& range)
{
    return static_cast<std::size_t>(std::distance(range.begin().base(), range.end().base()));
}

template <typename Iterator>
auto size(const ContainedRange<Iterator>& range, ForwardSortedTag)
{
    return static_cast<std::size_t>(std::distance(range.begin(), range.end()));
}

template <typename Iterator>
auto size(const ContainedRange<Iterator>& range, BidirectionallySortedTag)
{
    return base_size(range);
}

template <typename Iterator>
auto size(const ContainedRange<Iterator>& range)
{
    return size(range, ForwardSortedTag {});
}

template <typename Iterator>
bool empty(const ContainedRange<Iterator>& range)
{
    return range.empty();
}

template <typename Iterator, typename MappableType>
ContainedRange<Iterator> make_contained_range(Iterator first, Iterator last, const MappableType& mappable)
{
    using boost::make_iterator_range; using boost::make_filter_iterator; using detail::IsContained;
    using MappableType2 = typename std::iterator_traits<Iterator>::value_type;
    return make_iterator_range(
        make_filter_iterator<IsContained<MappableType2>>(IsContained<MappableType2>(mappable), first, last),
        make_filter_iterator<IsContained<MappableType2>>(IsContained<MappableType2>(mappable), last, last)
    );
}

namespace detail {

template <typename MappableType>
class IsShared
{
public:
    IsShared() = delete;
    
    template <typename MappableType1_, typename MappableType2_>
    IsShared(const MappableType1_& lhs, MappableType2_ rhs)
    : lhs_ {mapped_region(lhs)}
    , rhs_ {mapped_region(rhs)}
    {}
    
    bool operator()(const MappableType& mappable) const
    {
        return overlaps(lhs_, mappable) && overlaps(mappable, rhs_);
    }
    
private:
    RegionType<MappableType> lhs_, rhs_;
};

} // namespace detail

template <typename Iterator>
using SharedIterator = boost::filter_iterator<detail::IsShared<
    typename std::iterator_traits<Iterator>::value_type>,
    Iterator
>;

template <typename Iterator>
using SharedRange = boost::iterator_range<SharedIterator<Iterator>>;

template <typename Iterator>
boost::iterator_range<Iterator> bases(const SharedRange<Iterator>& range)
{
    return boost::make_iterator_range(range.begin().base(), range.end().base());
}

template <typename Iterator>
auto base_size(const SharedRange<Iterator>& range)
{
    return static_cast<std::size_t>(std::distance(range.begin().base(), range.end().base()));
}

template <typename Iterator>
auto size(const SharedRange<Iterator>& range, ForwardSortedTag)
{
    return static_cast<std::size_t>(std::distance(range.begin(), range.end()));
}

template <typename Iterator>
auto size(const SharedRange<Iterator>& range, BidirectionallySortedTag)
{
    return base_size(range);
}

template <typename Iterator>
auto size(const SharedRange<Iterator>& range)
{
    return size(range, ForwardSortedTag {});
}

template <typename Iterator>
bool empty(const SharedRange<Iterator>& range)
{
    return range.empty();
}

template <typename Iterator, typename MappableType1, typename MappableType2>
SharedRange<Iterator> make_shared_range(Iterator first, Iterator last,
                                        const MappableType1& mappable1, const MappableType2& mappable2)
{
    using boost::make_iterator_range; using boost::make_filter_iterator; using detail::IsShared;
    using MappableType_ = typename std::iterator_traits<Iterator>::value_type;
    return make_iterator_range(
        make_filter_iterator<IsShared<MappableType_>>(IsShared<MappableType_>(mappable1, mappable2), first, last),
        make_filter_iterator<IsShared<MappableType_>>(IsShared<MappableType_>(mappable1, mappable2), last, last)
    );
}

template <typename T>
auto crbegin(const boost::iterator_range<T>& range) noexcept
{
    return boost::rbegin(range);
}
    
template <typename T>
auto crend(const boost::iterator_range<T>& range) noexcept
{
    return boost::rend(range);
}

} // namespace octopus

#endif
