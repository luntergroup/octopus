//
//  mappable_set.h
//  Octopus
//
//  Created by Daniel Cooke on 19/07/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_mappable_set_h
#define Octopus_mappable_set_h

#include <memory>     // std::allocator
#include <functional> // std::less
#include <algorithm>  // std::max
#include <iterator>   // std::begin, std::end, std::cbegin, std::cend
#include <boost/container/flat_set.hpp>

#include "genomic_region.h"
#include "mappable.h"
#include "region_algorithms.h"

template <typename MappableType, typename Allocator = std::allocator<MappableType>>
class MappableSet
{
protected:
    using base_t = boost::container::flat_multiset<MappableType, std::less<MappableType>, Allocator>;
    
public:
    using allocator_type  = typename base_t::allocator_type;
    using value_type      = typename base_t::value_type;
    using reference       = typename base_t::reference ;
    using const_reference = typename base_t::const_reference;
    using difference_type = typename base_t::difference_type ;
    using size_type       = typename base_t::size_type ;
    
    using iterator               = typename base_t::iterator;
    using const_iterator         = typename base_t::const_iterator;
    using reverse_iterator       = typename base_t::reverse_iterator;
    using const_reverse_iterator = typename base_t::const_reverse_iterator;
    
    MappableSet() = default;
    template <typename InputIterator>
    MappableSet(InputIterator first, InputIterator second);
    ~MappableSet() = default;
    
    MappableSet(const MappableSet&)            = default;
    MappableSet& operator=(const MappableSet&) = default;
    MappableSet(MappableSet&&)                 = default;
    MappableSet& operator=(MappableSet&&)      = default;
    
    iterator begin();
    const_iterator begin() const;
    const_iterator cbegin() const;
    iterator end();
    const_iterator end() const;
    const_iterator cend() const;
    reverse_iterator rbegin();
    const_reverse_iterator rbegin() const;
    const_reverse_iterator crbegin() const;
    reverse_iterator rend();
    const_reverse_iterator rend() const;
    const_reverse_iterator crend() const;
    
    template <typename ...Args>
    iterator emplace(Args...);
    iterator insert(const MappableType&);
    iterator insert(MappableType&&);
    template <typename InputIterator>
    void insert(InputIterator, InputIterator);
    iterator insert(std::initializer_list<MappableType>);
    iterator erase(const_iterator);
    iterator erase(const MappableType&);
    iterator erase(const_iterator, const_iterator);
    void clear();
    
    void swap(const MappableSet&);
    size_type size() const noexcept;
    size_type capacity() const noexcept;
    size_type max_size() const noexcept;
    bool empty() const noexcept;
    void reserve(size_type n);
    void shrink_to_fit();
    
    allocator_type get_allocator();
    
    template <typename MappableType_>
    bool has_overlapped(const MappableType_& mappable);
    template <typename MappableType_>
    bool has_overlapped(iterator first, iterator last, const MappableType_& mappable);
    template <typename MappableType_>
    bool has_overlapped(const_iterator first, const_iterator last, const MappableType_& mappable);
    
    template <typename MappableType_>
    size_type count_overlapped(const MappableType_& mappable);
    template <typename MappableType_>
    size_type count_overlapped(iterator first, iterator last, const MappableType_& mappable);
    template <typename MappableType_>
    size_type count_overlapped(const_iterator first, const_iterator last, const MappableType_& mappable);
    
    template <typename MappableType_>
    OverlapRange<iterator> overlap_range(const MappableType_& mappable);
    template <typename MappableType_>
    OverlapRange<iterator> overlap_range(iterator first, iterator last, const MappableType_& mappable);
    template <typename MappableType_>
    OverlapRange<const_iterator> overlap_range(const_iterator first, const_iterator last, const MappableType_& mappable);
    
    template <typename MappableType_>
    bool has_contained(const MappableType_& mappable);
    template <typename MappableType_>
    bool has_contained(iterator first, iterator last, const MappableType_& mappable);
    template <typename MappableType_>
    bool has_contained(const_iterator first, const_iterator last, const MappableType_& mappable);
    
    template <typename MappableType_>
    size_type count_contained(const MappableType_& mappable);
    template <typename MappableType_>
    size_type count_contained(iterator first, iterator last, const MappableType_& mappable);
    template <typename MappableType_>
    size_type count_contained(const_iterator first, const_iterator last, const MappableType_& mappable);
    
    template <typename MappableType_>
    ContainedRange<iterator> contained_range(const MappableType_& mappable);
    template <typename MappableType_>
    ContainedRange<iterator> contained_range(iterator first, iterator last, const MappableType_& mappable);
    template <typename MappableType_>
    ContainedRange<const_iterator> contained_range(const_iterator first, const_iterator last, const MappableType_& mappable);
    
private:
    base_t elements_;
    
    bool is_bidirectionally_sorted_;
    GenomicRegion::SizeType max_element_size_;
};

template <typename MappableType, typename Allocator>
template <typename InputIterator>
MappableSet<MappableType, Allocator>::MappableSet(InputIterator first, InputIterator second)
:
elements_ {first, second},
is_bidirectionally_sorted_ {is_bidirectionally_sorted(std::cbegin(elements_), std::cend(elements_))},
max_element_size_ {::size(*largest(std::cbegin(elements_), std::cend(elements_)))}
{}

template <typename MappableType, typename Allocator>
typename MappableSet<MappableType, Allocator>::iterator
MappableSet<MappableType, Allocator>::begin()
{
    return elements_.begin();
}

template <typename MappableType, typename Allocator>
typename MappableSet<MappableType, Allocator>::const_iterator
MappableSet<MappableType, Allocator>::begin() const
{
    return elements_.begin();
}

template <typename MappableType, typename Allocator>
typename MappableSet<MappableType, Allocator>::const_iterator
MappableSet<MappableType, Allocator>::cbegin() const
{
    return elements_.cbegin();
}

template <typename MappableType, typename Allocator>
typename MappableSet<MappableType, Allocator>::iterator
MappableSet<MappableType, Allocator>::end()
{
    return elements_.end();
}

template <typename MappableType, typename Allocator>
typename MappableSet<MappableType, Allocator>::const_iterator
MappableSet<MappableType, Allocator>::end() const
{
    return elements_.end();
}

template <typename MappableType, typename Allocator>
typename MappableSet<MappableType, Allocator>::const_iterator
MappableSet<MappableType, Allocator>::cend() const
{
    return elements_.cend();
}

template <typename MappableType, typename Allocator>
typename MappableSet<MappableType, Allocator>::reverse_iterator
MappableSet<MappableType, Allocator>::rbegin()
{
    return elements_.rbegin();
}

template <typename MappableType, typename Allocator>
typename MappableSet<MappableType, Allocator>::const_reverse_iterator
MappableSet<MappableType, Allocator>::rbegin() const
{
    return elements_.rbegin();
}

template <typename MappableType, typename Allocator>
typename MappableSet<MappableType, Allocator>::const_reverse_iterator
MappableSet<MappableType, Allocator>::crbegin() const
{
    return elements_.crbegin();
}

template <typename MappableType, typename Allocator>
typename MappableSet<MappableType, Allocator>::reverse_iterator
MappableSet<MappableType, Allocator>::rend()
{
    return elements_.rend();
}

template <typename MappableType, typename Allocator>
typename MappableSet<MappableType, Allocator>::const_reverse_iterator
MappableSet<MappableType, Allocator>::rend() const
{
    return elements_.rend();
}

template <typename MappableType, typename Allocator>
typename MappableSet<MappableType, Allocator>::const_reverse_iterator
MappableSet<MappableType, Allocator>::crend() const
{
    return elements_.crend();
}

template <typename MappableType, typename Allocator>
template <typename ...Args>
typename MappableSet<MappableType, Allocator>::iterator
MappableSet<MappableType, Allocator>::emplace(Args... args)
{
    auto it = elements_.emplace(std::forward<Args>(args)...);
    
    if (it == std::prev(std::end(elements_))) {
        
    }
    
    max_element_size_ = std::max(max_element_size_, ::size(*it));
    return it;
}

template <typename MappableType, typename Allocator>
typename MappableSet<MappableType, Allocator>::iterator
MappableSet<MappableType, Allocator>::insert(const MappableType& m)
{
    auto it = elements_.insert(m);
    max_element_size_ = std::max(max_element_size_, ::size(*it));
    return it;
}

template <typename MappableType, typename Allocator>
typename MappableSet<MappableType, Allocator>::iterator
MappableSet<MappableType, Allocator>::insert(MappableType&& m)
{
    auto it = elements_.insert(std::move(m));
    max_element_size_ = std::max(max_element_size_, ::size(*it));
    return it;
}

template <typename MappableType, typename Allocator>
template <typename InputIterator>
void
MappableSet<MappableType, Allocator>::insert(InputIterator first, InputIterator last)
{
    max_element_size_ = std::max(max_element_size_, ::size(*largest(first, last)));
    elements_.insert(first, last);
}

template <typename MappableType, typename Allocator>
typename MappableSet<MappableType, Allocator>::iterator
MappableSet<MappableType, Allocator>::insert(std::initializer_list<MappableType> il)
{
    max_element_size_ = std::max(max_element_size_, ::size(*largest(il.begin(), il.end())));
    return elements_.insert(std::move(il));
}

template <typename MappableType, typename Allocator>
typename MappableSet<MappableType, Allocator>::iterator
MappableSet<MappableType, Allocator>::erase(const_iterator p)
{
    if (max_element_size_ == ::size(*p)) {
        auto it = elements_.erase(p);
        max_element_size_ = ::size(*largest(std::cbegin(elements_), std::cend(elements_)));
        return it;
    }
    return elements_.erase(p);
}

template <typename MappableType, typename Allocator>
typename MappableSet<MappableType, Allocator>::iterator
MappableSet<MappableType, Allocator>::erase(const MappableType& m)
{
    if (max_element_size_ == ::size(m)) {
        auto it = elements_.erase(m);
        max_element_size_ = ::size(*largest(std::cbegin(elements_), std::cend(elements_)));
        return it;
    }
    return elements_.erase(m);
}

template <typename MappableType, typename Allocator>
typename MappableSet<MappableType, Allocator>::iterator
MappableSet<MappableType, Allocator>::erase(const_iterator first, const_iterator last)
{
    if (max_element_size_ == ::size(*largest(first, last))) {
        auto it = elements_.erase(first, last);
        max_element_size_ = ::size(*largest(std::cbegin(elements_), std::cend(elements_)));
        return it;
    }
    return elements_.erase(first, last);
}

template <typename MappableType, typename Allocator>
void MappableSet<MappableType, Allocator>::clear()
{
    elements_.clear();
    max_element_size_ = 0;
    is_bidirectionally_sorted_ = true;
}

template <typename MappableType, typename Allocator>
void MappableSet<MappableType, Allocator>::swap(const MappableSet& m)
{
    elements_.swap(m.elements_);
}

template <typename MappableType, typename Allocator>
typename MappableSet<MappableType, Allocator>::size_type
MappableSet<MappableType, Allocator>::size() const noexcept
{
    return elements_.size();
}

template <typename MappableType, typename Allocator>
typename MappableSet<MappableType, Allocator>::size_type
MappableSet<MappableType, Allocator>::capacity() const noexcept
{
    return elements_.capacity();
}

template <typename MappableType, typename Allocator>
typename MappableSet<MappableType, Allocator>::size_type
MappableSet<MappableType, Allocator>::max_size() const noexcept
{
    return elements_.max_size();
}

template <typename MappableType, typename Allocator>
bool MappableSet<MappableType, Allocator>::empty() const noexcept
{
    return elements_.empty();
}

template <typename MappableType, typename Allocator>
void
MappableSet<MappableType, Allocator>::reserve(size_type n)
{
    elements_.reserve(n);
}

template <typename MappableType, typename Allocator>
void
MappableSet<MappableType, Allocator>::shrink_to_fit()
{
    elements_.shrink_to_fit();
}

template <typename MappableType, typename Allocator>
typename MappableSet<MappableType, Allocator>::allocator_type
MappableSet<MappableType, Allocator>::get_allocator()
{
    return elements_.get_allocator();
}

template <typename MappableType, typename Allocator>
template <typename MappableType_>
bool
MappableSet<MappableType, Allocator>::has_overlapped(const MappableType_& mappable)
{
    return (is_bidirectionally_sorted_) ?
    ::has_overlapped(std::begin(elements_), std::end(elements_), mappable, MappableRangeOrder::BidirectionallySorted)
    :
    ::has_overlapped(std::begin(elements_), std::end(elements_), mappable);
}

template <typename MappableType, typename Allocator>
template <typename MappableType_>
bool
MappableSet<MappableType, Allocator>::has_overlapped(iterator first, iterator last, const MappableType_& mappable)
{
    return (is_bidirectionally_sorted_) ?
    ::has_overlapped(first, last, mappable, MappableRangeOrder::BidirectionallySorted)
    :
    ::has_overlapped(first, last, mappable, max_element_size_);
}

template <typename MappableType, typename Allocator>
template <typename MappableType_>
bool
MappableSet<MappableType, Allocator>::has_overlapped(const_iterator first, const_iterator last, const MappableType_& mappable)
{
    return (is_bidirectionally_sorted_) ?
    ::has_overlapped(first, last, mappable, MappableRangeOrder::BidirectionallySorted)
    :
    ::has_overlapped(first, last, mappable, max_element_size_);
}

template <typename MappableType, typename Allocator>
template <typename MappableType_>
typename MappableSet<MappableType, Allocator>::size_type
MappableSet<MappableType, Allocator>::count_overlapped(const MappableType_& mappable)
{
    auto overlapped = overlap_range(mappable);
    return (is_bidirectionally_sorted_) ? ::size(bases(overlapped)) : ::size(overlapped);
}

template <typename MappableType, typename Allocator>
template <typename MappableType_>
typename MappableSet<MappableType, Allocator>::size_type
MappableSet<MappableType, Allocator>::count_overlapped(iterator first, iterator last, const MappableType_& mappable)
{
    auto overlapped = overlap_range(first, last, mappable);
    return (is_bidirectionally_sorted_) ? ::size(bases(overlapped)) : ::size(overlapped);
}

template <typename MappableType, typename Allocator>
template <typename MappableType_>
typename MappableSet<MappableType, Allocator>::size_type
MappableSet<MappableType, Allocator>::count_overlapped(const_iterator first, const_iterator last, const MappableType_& mappable)
{
    auto overlapped = overlap_range(first, last, mappable);
    return (is_bidirectionally_sorted_) ? ::size(bases(overlapped)) : ::size(overlapped);
}

template <typename MappableType, typename Allocator>
template <typename MappableType_>
OverlapRange<typename MappableSet<MappableType, Allocator>::iterator>
MappableSet<MappableType, Allocator>::overlap_range(const MappableType_& mappable)
{
    return (is_bidirectionally_sorted_) ?
        ::overlap_range(std::begin(elements_), std::end(elements_), mappable, MappableRangeOrder::BidirectionallySorted)
    :
        ::overlap_range(std::begin(elements_), std::end(elements_), mappable, max_element_size_);
}

template <typename MappableType, typename Allocator>
template <typename MappableType_>
OverlapRange<typename MappableSet<MappableType, Allocator>::iterator>
MappableSet<MappableType, Allocator>::overlap_range(iterator first, iterator last, const MappableType_& mappable)
{
    return (is_bidirectionally_sorted_) ?
        ::overlap_range(first, last, mappable, MappableRangeOrder::BidirectionallySorted)
    :
        ::overlap_range(first, last, mappable, max_element_size_);
}

template <typename MappableType, typename Allocator>
template <typename MappableType_>
OverlapRange<typename MappableSet<MappableType, Allocator>::const_iterator>
MappableSet<MappableType, Allocator>::overlap_range(const_iterator first, const_iterator last, const MappableType_& mappable)
{
    return (is_bidirectionally_sorted_) ?
        ::overlap_range(first, last, mappable, MappableRangeOrder::BidirectionallySorted)
    :
        ::overlap_range(first, last, mappable, max_element_size_);
}

template <typename MappableType, typename Allocator>
template <typename MappableType_>
bool
MappableSet<MappableType, Allocator>::has_contained(const MappableType_& mappable)
{
    return ::has_contained(std::begin(elements_), std::end(elements_), mappable);
    
}

template <typename MappableType, typename Allocator>
template <typename MappableType_>
bool
MappableSet<MappableType, Allocator>::has_contained(iterator first, iterator last, const MappableType_& mappable)
{
    return ::has_contained(first, last, mappable);
}

template <typename MappableType, typename Allocator>
template <typename MappableType_>
bool
MappableSet<MappableType, Allocator>::has_contained(const_iterator first, const_iterator last, const MappableType_& mappable)
{
    return ::has_contained(first, last, mappable);
}

template <typename MappableType, typename Allocator>
template <typename MappableType_>
typename MappableSet<MappableType, Allocator>::size_type
MappableSet<MappableType, Allocator>::count_contained(const MappableType_& mappable)
{
    auto contained = contained_range(mappable);
    return (is_bidirectionally_sorted_) ? ::size(bases(contained)) : ::size(contained);
}

template <typename MappableType, typename Allocator>
template <typename MappableType_>
typename MappableSet<MappableType, Allocator>::size_type
MappableSet<MappableType, Allocator>::count_contained(iterator first, iterator last, const MappableType_& mappable)
{
    auto contained = contained_range(first, last, mappable);
    return (is_bidirectionally_sorted_) ? ::size(bases(contained)) : ::size(contained);
}

template <typename MappableType, typename Allocator>
template <typename MappableType_>
typename MappableSet<MappableType, Allocator>::size_type
MappableSet<MappableType, Allocator>::count_contained(const_iterator first, const_iterator last, const MappableType_& mappable)
{
    auto contained = contained_range(first, last, mappable);
    return (is_bidirectionally_sorted_) ? ::size(bases(contained)) : ::size(contained);
}

template <typename MappableType, typename Allocator>
template <typename MappableType_>
ContainedRange<typename MappableSet<MappableType, Allocator>::iterator>
MappableSet<MappableType, Allocator>::contained_range(const MappableType_& mappable)
{
    return ::contained_range(std::begin(elements_), std::end(elements_), mappable);
}

template <typename MappableType, typename Allocator>
template <typename MappableType_>
ContainedRange<typename MappableSet<MappableType, Allocator>::iterator>
MappableSet<MappableType, Allocator>::contained_range(iterator first, iterator last, const MappableType_& mappable)
{
    return ::contained_range(first, last, mappable);
}

template <typename MappableType, typename Allocator>
template <typename MappableType_>
ContainedRange<typename MappableSet<MappableType, Allocator>::const_iterator>
MappableSet<MappableType, Allocator>::contained_range(const_iterator first, const_iterator last, const MappableType_& mappable)
{
    return ::contained_range(first, last, mappable);
}

#endif
