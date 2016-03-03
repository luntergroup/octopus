//
//  mappable_flat_set.hpp
//  Octopus
//
//  Created by Daniel Cooke on 02/03/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#ifndef mappable_flat_set_hpp
#define mappable_flat_set_hpp

#include <deque>
#include <memory>
#include <functional>
#include <algorithm>
#include <iterator>
#include <initializer_list>
#include <stdexcept>
#include <vector>

#include "comparable.hpp"
#include "mappable.hpp"
#include "mappable_ranges.hpp"
#include "mappable_algorithms.hpp"

#include <iostream> // TEST

/*
 MappableFlatSet is a container designed to allow fast retrival of MappableType elements with minimal
 memory overhead.
 */
template <typename MappableType, typename Allocator = std::allocator<MappableType>>
class MappableFlatSet : public Comparable<MappableFlatSet<MappableType, Allocator>>
{
protected:
    using base_t = std::deque<MappableType, Allocator>;
    
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
    
    MappableFlatSet();
    
    template <typename InputIterator>
    MappableFlatSet(InputIterator first, InputIterator second);
    
    ~MappableFlatSet() = default;
    
    MappableFlatSet(const MappableFlatSet&)            = default;
    MappableFlatSet& operator=(const MappableFlatSet&) = default;
    MappableFlatSet(MappableFlatSet&&)                 = default;
    MappableFlatSet& operator=(MappableFlatSet&&)      = default;
    
    iterator begin() noexcept;
    const_iterator begin() const noexcept;
    const_iterator cbegin() const noexcept;
    iterator end() noexcept;
    const_iterator end() const noexcept;
    const_iterator cend() const noexcept;
    reverse_iterator rbegin() noexcept;
    const_reverse_iterator rbegin() const noexcept;
    const_reverse_iterator crbegin() const noexcept;
    reverse_iterator rend() noexcept;
    const_reverse_iterator rend() const noexcept;
    const_reverse_iterator crend() const noexcept;
    
    reference at(size_type pos);
    const_reference at(size_type pos) const;
    reference operator[](size_type pos);
    const_reference operator[](size_type pos) const;
    reference front();
    const_reference front() const;
    reference back();
    const_reference back() const;
    
//    template <typename ...Args>
//    iterator emplace(Args...);
//    iterator insert(const MappableType&);
//    iterator insert(MappableType&&);
//    iterator insert(const_iterator, const MappableType& mappable);
//    iterator insert(const_iterator, MappableType&& mappable);
//    template <typename InputIterator>
//    void insert(InputIterator, InputIterator);
//    iterator insert(std::initializer_list<MappableType>);
    iterator erase(const_iterator);
    size_type erase(const MappableType&);
    iterator erase(const_iterator, const_iterator);
    template <typename InputIt>
    size_type erase_all(InputIt first, InputIt last);
    
    void clear();
    
    void swap(const MappableFlatSet&);
    
    size_type size() const noexcept;
    size_type capacity() const noexcept;
    size_type max_size() const noexcept;
    bool empty() const noexcept;
    void shrink_to_fit();
    
    const MappableType& leftmost() const;
    const MappableType& rightmost() const;
    
    template <typename MappableType_>
    bool has_overlapped(const MappableType_& mappable) const;
    template <typename MappableType_>
    bool has_overlapped(iterator first, iterator last, const MappableType_& mappable) const;
    template <typename MappableType_>
    bool has_overlapped(const_iterator first, const_iterator last, const MappableType_& mappable) const;
    
    template <typename MappableType_>
    size_type count_overlapped(const MappableType_& mappable) const;
    template <typename MappableType_>
    size_type count_overlapped(iterator first, iterator last, const MappableType_& mappable) const;
    template <typename MappableType_>
    size_type count_overlapped(const_iterator first, const_iterator last,
                               const MappableType_& mappable) const;
    
    template <typename MappableType_>
    OverlapRange<const_iterator> overlap_range(const MappableType_& mappable) const;
    template <typename MappableType_>
    OverlapRange<iterator> overlap_range(iterator first, iterator last,
                                         const MappableType_& mappable) const;
    template <typename MappableType_>
    OverlapRange<const_iterator> overlap_range(const_iterator first, const_iterator last,
                                               const MappableType_& mappable) const;
    
    template <typename MappableType_>
    void erase_overlapped(const MappableType_& mappable);
    
    template <typename MappableType_>
    bool has_contained(const MappableType_& mappable) const;
    template <typename MappableType_>
    bool has_contained(iterator first, iterator last, const MappableType_& mappable) const;
    template <typename MappableType_>
    bool has_contained(const_iterator first, const_iterator last, const MappableType_& mappable) const;
    
    template <typename MappableType_>
    size_type count_contained(const MappableType_& mappable) const;
    template <typename MappableType_>
    size_type count_contained(iterator first, iterator last, const MappableType_& mappable) const;
    template <typename MappableType_>
    size_type count_contained(const_iterator first, const_iterator last,
                              const MappableType_& mappable) const;
    
    template <typename MappableType_>
    ContainedRange<const_iterator> contained_range(const MappableType_& mappable) const;
    template <typename MappableType_>
    ContainedRange<iterator> contained_range(iterator first, iterator last,
                                             const MappableType_& mappable) const;
    template <typename MappableType_>
    ContainedRange<const_iterator> contained_range(const_iterator first, const_iterator last,
                                                   const MappableType_& mappable) const;
    
    template <typename MappableType_>
    void erase_contained(const MappableType_& mappable);
    
    template <typename M, typename A>
    friend bool operator==(const MappableFlatSet<M, A>& lhs, const MappableFlatSet<M, A>& rhs);
    template <typename M, typename A>
    friend bool operator<(const MappableFlatSet<M, A>& lhs, const MappableFlatSet<M, A>& rhs);
    template <typename M, typename A>
    friend void swap(MappableFlatSet<M, A>& lhs, MappableFlatSet<M, A>& rhs);
    
private:
    base_t elements_;
    bool is_bidirectionally_sorted_;
    typename RegionType<MappableType>::SizeType max_element_size_;
};

template <typename MappableType, typename Allocator>
MappableFlatSet<MappableType, Allocator>::MappableFlatSet()
:
elements_ {},
is_bidirectionally_sorted_ {true},
max_element_size_ {}
{}

template <typename MappableType, typename Allocator>
template <typename InputIterator>
MappableFlatSet<MappableType, Allocator>::MappableFlatSet(InputIterator first, InputIterator second)
:
elements_ {first, second},
is_bidirectionally_sorted_ {true},
max_element_size_ {(elements_.empty()) ? 0 : region_size(*largest_mappable(elements_))}
{
    std::sort(std::begin(elements_), std::end(elements_));
    
    elements_.erase(std::unique(std::begin(elements_), std::end(elements_)), std::end(elements_));
    
    is_bidirectionally_sorted_ = is_bidirectionally_sorted(elements_);
}

template <typename MappableType, typename Allocator>
typename MappableFlatSet<MappableType, Allocator>::iterator
MappableFlatSet<MappableType, Allocator>::begin() noexcept
{
    return elements_.begin();
}

template <typename MappableType, typename Allocator>
typename MappableFlatSet<MappableType, Allocator>::const_iterator
MappableFlatSet<MappableType, Allocator>::begin() const noexcept
{
    return elements_.begin();
}

template <typename MappableType, typename Allocator>
typename MappableFlatSet<MappableType, Allocator>::const_iterator
MappableFlatSet<MappableType, Allocator>::cbegin() const noexcept
{
    return elements_.cbegin();
}

template <typename MappableType, typename Allocator>
typename MappableFlatSet<MappableType, Allocator>::iterator
MappableFlatSet<MappableType, Allocator>::end() noexcept
{
    return elements_.end();
}

template <typename MappableType, typename Allocator>
typename MappableFlatSet<MappableType, Allocator>::const_iterator
MappableFlatSet<MappableType, Allocator>::end() const noexcept
{
    return elements_.end();
}

template <typename MappableType, typename Allocator>
typename MappableFlatSet<MappableType, Allocator>::const_iterator
MappableFlatSet<MappableType, Allocator>::cend() const noexcept
{
    return elements_.cend();
}

template <typename MappableType, typename Allocator>
typename MappableFlatSet<MappableType, Allocator>::reverse_iterator
MappableFlatSet<MappableType, Allocator>::rbegin() noexcept
{
    return elements_.rbegin();
}

template <typename MappableType, typename Allocator>
typename MappableFlatSet<MappableType, Allocator>::const_reverse_iterator
MappableFlatSet<MappableType, Allocator>::rbegin() const noexcept
{
    return elements_.rbegin();
}

template <typename MappableType, typename Allocator>
typename MappableFlatSet<MappableType, Allocator>::const_reverse_iterator
MappableFlatSet<MappableType, Allocator>::crbegin() const noexcept
{
    return elements_.crbegin();
}

template <typename MappableType, typename Allocator>
typename MappableFlatSet<MappableType, Allocator>::reverse_iterator
MappableFlatSet<MappableType, Allocator>::rend() noexcept
{
    return elements_.rend();
}

template <typename MappableType, typename Allocator>
typename MappableFlatSet<MappableType, Allocator>::const_reverse_iterator
MappableFlatSet<MappableType, Allocator>::rend() const noexcept
{
    return elements_.rend();
}

template <typename MappableType, typename Allocator>
typename MappableFlatSet<MappableType, Allocator>::const_reverse_iterator
MappableFlatSet<MappableType, Allocator>::crend() const noexcept
{
    return elements_.crend();
}

template <typename MappableType, typename Allocator>
typename MappableFlatSet<MappableType, Allocator>::reference
MappableFlatSet<MappableType, Allocator>::at(size_type pos)
{
    if (pos < size()) {
        return *std::next(begin(), pos);
    } else {
        throw std::out_of_range {"MappableFlatSet"};
    }
}

template <typename MappableType, typename Allocator>
typename MappableFlatSet<MappableType, Allocator>::const_reference
MappableFlatSet<MappableType, Allocator>::at(size_type pos) const
{
    if (pos < size()) {
        return *std::next(cbegin(), pos);
    } else {
        throw std::out_of_range {"MappableFlatSet"};
    }
}

template <typename MappableType, typename Allocator>
typename MappableFlatSet<MappableType, Allocator>::reference
MappableFlatSet<MappableType, Allocator>::operator[](size_type pos)
{
    return *std::next(begin(), pos);
}

template <typename MappableType, typename Allocator>
typename MappableFlatSet<MappableType, Allocator>::const_reference
MappableFlatSet<MappableType, Allocator>::operator[](size_type pos) const
{
    return *std::next(cbegin(), pos);
}

template <typename MappableType, typename Allocator>
typename MappableFlatSet<MappableType, Allocator>::reference
MappableFlatSet<MappableType, Allocator>::front()
{
    return *begin();
}

template <typename MappableType, typename Allocator>
typename MappableFlatSet<MappableType, Allocator>::const_reference
MappableFlatSet<MappableType, Allocator>::front() const
{
    return *cbegin();
}

template <typename MappableType, typename Allocator>
typename MappableFlatSet<MappableType, Allocator>::reference
MappableFlatSet<MappableType, Allocator>::back()
{
    return *std::prev(end());
}

template <typename MappableType, typename Allocator>
typename MappableFlatSet<MappableType, Allocator>::const_reference
MappableFlatSet<MappableType, Allocator>::back() const
{
    return *std::prev(cend());
}

//template <typename MappableType, typename Allocator>
//template <typename ...Args>
//typename MappableFlatSet<MappableType, Allocator>::iterator
//MappableFlatSet<MappableType, Allocator>::emplace(Args... args)
//{
//    const auto it = elements_.emplace(std::forward<Args>(args)...);
//    if (is_bidirectionally_sorted_) {
//        const auto overlapped = overlap_range(*it);
//        is_bidirectionally_sorted_ = is_bidirectionally_sorted(overlapped);
//    }
//    max_element_size_ = std::max(max_element_size_, region_size(*it));
//    return it;
//}
//
//template <typename MappableType, typename Allocator>
//typename MappableFlatSet<MappableType, Allocator>::iterator
//MappableFlatSet<MappableType, Allocator>::insert(const MappableType& m)
//{
//    const auto it = elements_.insert(m);
//    if (is_bidirectionally_sorted_) {
//        const auto overlapped = overlap_range(*it);
//        is_bidirectionally_sorted_ = is_bidirectionally_sorted(overlapped);
//    }
//    max_element_size_ = std::max(max_element_size_, region_size(*it));
//    return it;
//}
//
//template <typename MappableType, typename Allocator>
//typename MappableFlatSet<MappableType, Allocator>::iterator
//MappableFlatSet<MappableType, Allocator>::insert(MappableType&& m)
//{
//    const auto it = elements_.insert(std::move(m));
//    if (is_bidirectionally_sorted_) {
//        const auto overlapped = overlap_range(*it);
//        is_bidirectionally_sorted_ = is_bidirectionally_sorted(overlapped);
//    }
//    max_element_size_ = std::max(max_element_size_, region_size(*it));
//    return it;
//}
//
//template <typename MappableType, typename Allocator>
//typename MappableFlatSet<MappableType, Allocator>::iterator
//MappableFlatSet<MappableType, Allocator>::insert(const_iterator it, const MappableType& m)
//{
//    const auto it2 = elements_.insert(it, m);
//    if (is_bidirectionally_sorted_) {
//        const auto overlapped = overlap_range(*it2);
//        is_bidirectionally_sorted_ = is_bidirectionally_sorted(overlapped);
//    }
//    max_element_size_ = std::max(max_element_size_, region_size(*it));
//    return it2;
//}
//
//template <typename MappableType, typename Allocator>
//typename MappableFlatSet<MappableType, Allocator>::iterator
//MappableFlatSet<MappableType, Allocator>::insert(const_iterator it, MappableType&& m)
//{
//    const auto it2 = elements_.insert(it, std::move(m));
//    if (is_bidirectionally_sorted_) {
//        const auto overlapped = overlap_range(*it2);
//        is_bidirectionally_sorted_ = is_bidirectionally_sorted(overlapped);
//    }
//    max_element_size_ = std::max(max_element_size_, region_size(*it));
//    return it2;
//}
//
//template <typename MappableType, typename Allocator>
//template <typename InputIterator>
//void
//MappableFlatSet<MappableType, Allocator>::insert(InputIterator first, InputIterator last)
//{
//    if (first != last) {
//        max_element_size_ = std::max(max_element_size_, region_size(*largest_mappable(first, last)));
//    }
//    elements_.insert(first, last);
//    if (is_bidirectionally_sorted_) {
//        is_bidirectionally_sorted_ = is_bidirectionally_sorted(elements_);
//    }
//}
//
//template <typename MappableType, typename Allocator>
//typename MappableFlatSet<MappableType, Allocator>::iterator
//MappableFlatSet<MappableType, Allocator>::insert(std::initializer_list<MappableType> il)
//{
//    if (!il.empty()) {
//        max_element_size_ = std::max(max_element_size_, region_size(*largest_element(il)));
//    }
//    return elements_.insert(std::move(il));
//    if (is_bidirectionally_sorted_) {
//        is_bidirectionally_sorted_ = is_bidirectionally_sorted(elements_);
//    }
//}

template <typename MappableType, typename Allocator>
typename MappableFlatSet<MappableType, Allocator>::iterator
MappableFlatSet<MappableType, Allocator>::erase(const_iterator p)
{
    if (p == cend()) return elements_.erase(p);
    
    const auto erased_size = region_size(*p);
    
    const auto result = elements_.erase(p);
    
//    if (!is_bidirectionally_sorted_) {
//        is_bidirectionally_sorted_ = is_bidirectionally_sorted(elements_);
//    }
    
    if (max_element_size_ == erased_size) {
        max_element_size_ = region_size(*largest_mappable(elements_));
    }
    
    return result;
}

template <typename MappableType, typename Allocator>
typename MappableFlatSet<MappableType, Allocator>::size_type
MappableFlatSet<MappableType, Allocator>::erase(const MappableType& m)
{
    const auto contained = contained_range(m);
    
    const auto it = std::find(std::cbegin(contained), std::cend(contained), m);
    
    if (it != std::cend(contained)) {
        erase(it.base());
        return 1;
    }
    
    return 0;
}

template <typename MappableType, typename Allocator>
typename MappableFlatSet<MappableType, Allocator>::iterator
MappableFlatSet<MappableType, Allocator>::erase(const_iterator first, const_iterator last)
{
    if (first == last) return elements_.erase(first, last);
    
    const auto max_erased_size = region_size(*largest_mappable(first, last));
    
    const auto result = elements_.erase(first, last);
    
//    if (!is_bidirectionally_sorted_) {
//        is_bidirectionally_sorted_ = is_bidirectionally_sorted(elements_);
//    }
    
    if (max_element_size_ == max_erased_size) {
        max_element_size_ = region_size(*largest_mappable(elements_));
    }
    
    return result;
}

template <typename MappableType, typename Allocator>
template <typename InputIt>
typename MappableFlatSet<MappableType, Allocator>::size_type
MappableFlatSet<MappableType, Allocator>::erase_all(InputIt first, InputIt last)
{
    size_type result {0};
    
    if (first == last) return result;
    
    const auto region = encompassing_region(first, last);
    
    const auto p = bases(contained_range(region));
    
    auto from     = std::cbegin(p);
    const auto to = std::cend(p);
    
    typename RegionType<MappableType>::SizeType max_erased_size {0};
    
    while (first != last) {
        const auto contained = contained_range(from, to, *first);
        
        const auto it = std::find(std::cbegin(contained), std::cend(contained), *first);
        
        if (it != std::cend(contained)) {
            from = erase(it.base());
            
            if (region_size(*first) > max_erased_size) {
                max_erased_size = region_size(*first);
            }
            
            ++result;
        }
        
        ++first;
    }
    
    if (result > 0) {
//        if (!is_bidirectionally_sorted_) {
//            is_bidirectionally_sorted_ = is_bidirectionally_sorted(elements_);
//        }
        
        if (max_element_size_ == max_erased_size) {
            max_element_size_ = region_size(*largest_mappable(elements_));
        }
    }
    
    return 0;
}

template <typename MappableType, typename Allocator>
void MappableFlatSet<MappableType, Allocator>::clear()
{
    elements_.clear();
    is_bidirectionally_sorted_ = true;
    max_element_size_ = 0;
}

template <typename MappableType, typename Allocator>
void MappableFlatSet<MappableType, Allocator>::swap(const MappableFlatSet& m)
{
    std::swap(elements_, m.elements_);
    std::swap(is_bidirectionally_sorted_, m.is_bidirectionally_sorted_);
    std::swap(max_element_size_, m.max_element_size_);
}

template <typename MappableType, typename Allocator>
typename MappableFlatSet<MappableType, Allocator>::size_type
MappableFlatSet<MappableType, Allocator>::size() const noexcept
{
    return elements_.size();
}

template <typename MappableType, typename Allocator>
typename MappableFlatSet<MappableType, Allocator>::size_type
MappableFlatSet<MappableType, Allocator>::max_size() const noexcept
{
    return elements_.max_size();
}

template <typename MappableType, typename Allocator>
bool MappableFlatSet<MappableType, Allocator>::empty() const noexcept
{
    return elements_.empty();
}

template <typename MappableType, typename Allocator>
void
MappableFlatSet<MappableType, Allocator>::shrink_to_fit()
{
    elements_.shrink_to_fit();
}

template <typename MappableType, typename Allocator>
const MappableType& MappableFlatSet<MappableType, Allocator>::leftmost() const
{
    return front();
}

template <typename MappableType, typename Allocator>
const MappableType& MappableFlatSet<MappableType, Allocator>::rightmost() const
{
    using std::cbegin; using std::cend;
    
    const auto& last = *std::prev(elements_.cend());
    
    if (is_bidirectionally_sorted_) {
        return last;
    } else {
        const auto overlapped = ::overlap_range(cbegin(elements_), cend(elements_), last,
                                                max_element_size_);
        return *rightmost_mappable(cbegin(overlapped), cend(overlapped));
    }
}

template <typename MappableType, typename Allocator>
template <typename MappableType_>
bool
MappableFlatSet<MappableType, Allocator>::has_overlapped(const MappableType_& mappable) const
{
    return (is_bidirectionally_sorted_) ?
    ::has_overlapped(std::begin(elements_), std::end(elements_), mappable,
                     MappableRangeOrder::BidirectionallySorted)
    :
    ::has_overlapped(std::begin(elements_), std::end(elements_), mappable);
}

template <typename MappableType, typename Allocator>
template <typename MappableType_>
bool
MappableFlatSet<MappableType, Allocator>::has_overlapped(iterator first, iterator last,
                                                              const MappableType_& mappable) const
{
    return (is_bidirectionally_sorted_) ?
    ::has_overlapped(first, last, mappable, MappableRangeOrder::BidirectionallySorted)
    :
    ::has_overlapped(first, last, mappable, max_element_size_);
}

template <typename MappableType, typename Allocator>
template <typename MappableType_>
bool
MappableFlatSet<MappableType, Allocator>::has_overlapped(const_iterator first, const_iterator last,
                                                              const MappableType_& mappable) const
{
    return (is_bidirectionally_sorted_) ?
    ::has_overlapped(first, last, mappable, MappableRangeOrder::BidirectionallySorted)
    :
    ::has_overlapped(first, last, mappable, max_element_size_);
}

template <typename MappableType, typename Allocator>
template <typename MappableType_>
typename MappableFlatSet<MappableType, Allocator>::size_type
MappableFlatSet<MappableType, Allocator>::count_overlapped(const MappableType_& mappable) const
{
    const auto overlapped = overlap_range(mappable);
    return (is_bidirectionally_sorted_) ? ::size(overlapped, MappableRangeOrder::BidirectionallySorted) : ::size(overlapped);
}

template <typename MappableType, typename Allocator>
template <typename MappableType_>
typename MappableFlatSet<MappableType, Allocator>::size_type
MappableFlatSet<MappableType, Allocator>::count_overlapped(iterator first, iterator last,
                                                                const MappableType_& mappable) const
{
    const auto overlapped = overlap_range(first, last, mappable);
    return (is_bidirectionally_sorted_) ? ::size(overlapped, MappableRangeOrder::BidirectionallySorted) : ::size(overlapped);
}

template <typename MappableType, typename Allocator>
template <typename MappableType_>
typename MappableFlatSet<MappableType, Allocator>::size_type
MappableFlatSet<MappableType, Allocator>::count_overlapped(const_iterator first, const_iterator last,
                                                                const MappableType_& mappable) const
{
    const auto overlapped = overlap_range(first, last, mappable);
    return (is_bidirectionally_sorted_) ? ::size(overlapped, MappableRangeOrder::BidirectionallySorted) : ::size(overlapped);
}

template <typename MappableType, typename Allocator>
template <typename MappableType_>
OverlapRange<typename MappableFlatSet<MappableType, Allocator>::const_iterator>
MappableFlatSet<MappableType, Allocator>::overlap_range(const MappableType_& mappable) const
{
    return overlap_range(std::cbegin(elements_), std::cend(elements_), mappable);
}

template <typename MappableType, typename Allocator>
template <typename MappableType_>
OverlapRange<typename MappableFlatSet<MappableType, Allocator>::iterator>
MappableFlatSet<MappableType, Allocator>::overlap_range(iterator first, iterator last,
                                                             const MappableType_& mappable) const
{
    return overlap_range(const_iterator(first), const_iterator(last), mappable);
}

template <typename MappableType, typename Allocator>
template <typename MappableType_>
OverlapRange<typename MappableFlatSet<MappableType, Allocator>::const_iterator>
MappableFlatSet<MappableType, Allocator>::overlap_range(const_iterator first, const_iterator last,
                                                             const MappableType_& mappable) const
{
    return (is_bidirectionally_sorted_) ?
    ::overlap_range(first, last, mappable, MappableRangeOrder::BidirectionallySorted)
    :
    ::overlap_range(first, last, mappable, max_element_size_);
}

template <typename MappableType, typename Allocator>
template <typename MappableType_>
void MappableFlatSet<MappableType, Allocator>::erase_overlapped(const MappableType_& mappable)
{
    // TODO: find better implementation
    
    const auto overlapped = overlap_range(mappable);
    
    if (is_bidirectionally_sorted_ || ::size(overlapped) == bases(overlapped).size()) {
        erase(std::cbegin(overlapped).base(), std::cend(overlapped).base());
    } else {
        const std::vector<MappableType> tmp {std::cbegin(overlapped), std::cend(overlapped)};
        
        erase_all(std::cbegin(tmp), std::cend(tmp));
    }
}

template <typename MappableType, typename Allocator>
template <typename MappableType_>
bool
MappableFlatSet<MappableType, Allocator>::has_contained(const MappableType_& mappable) const
{
    return has_contained(std::cbegin(elements_), std::cend(elements_), mappable);
}

template <typename MappableType, typename Allocator>
template <typename MappableType_>
bool
MappableFlatSet<MappableType, Allocator>::has_contained(iterator first, iterator last,
                                                             const MappableType_& mappable) const
{
    return has_contained(const_iterator(first), const_iterator(last), mappable);
}

template <typename MappableType, typename Allocator>
template <typename MappableType_>
bool
MappableFlatSet<MappableType, Allocator>::has_contained(const_iterator first, const_iterator last,
                                                             const MappableType_& mappable) const
{
    return ::has_contained(first, last, mappable);
}

template <typename MappableType, typename Allocator>
template <typename MappableType_>
typename MappableFlatSet<MappableType, Allocator>::size_type
MappableFlatSet<MappableType, Allocator>::count_contained(const MappableType_& mappable) const
{
    return count_contained(std::cbegin(elements_), std::cend(elements_), mappable);
}

template <typename MappableType, typename Allocator>
template <typename MappableType_>
typename MappableFlatSet<MappableType, Allocator>::size_type
MappableFlatSet<MappableType, Allocator>::count_contained(iterator first, iterator last,
                                                               const MappableType_& mappable) const
{
    return count_contained(const_iterator(first), const_iterator(last), mappable);
}

template <typename MappableType, typename Allocator>
template <typename MappableType_>
typename MappableFlatSet<MappableType, Allocator>::size_type
MappableFlatSet<MappableType, Allocator>::count_contained(const_iterator first, const_iterator last,
                                                               const MappableType_& mappable) const
{
    const auto contained = contained_range(first, last, mappable);
    return (is_bidirectionally_sorted_) ? ::size(contained, MappableRangeOrder::BidirectionallySorted) : ::size(contained);
}

template <typename MappableType, typename Allocator>
template <typename MappableType_>
ContainedRange<typename MappableFlatSet<MappableType, Allocator>::const_iterator>
MappableFlatSet<MappableType, Allocator>::contained_range(const MappableType_& mappable) const
{
    return contained_range(std::cbegin(elements_), std::cend(elements_), mappable);
}

template <typename MappableType, typename Allocator>
template <typename MappableType_>
ContainedRange<typename MappableFlatSet<MappableType, Allocator>::iterator>
MappableFlatSet<MappableType, Allocator>::contained_range(iterator first, iterator last,
                                                               const MappableType_& mappable) const
{
    return contained_range(const_iterator(first), const_iterator(last), mappable);
}

template <typename MappableType, typename Allocator>
template <typename MappableType_>
ContainedRange<typename MappableFlatSet<MappableType, Allocator>::const_iterator>
MappableFlatSet<MappableType, Allocator>::contained_range(const_iterator first, const_iterator last,
                                                               const MappableType_& mappable) const
{
    return ::contained_range(first, last, mappable);
}

template <typename MappableType, typename Allocator>
template <typename MappableType_>
void MappableFlatSet<MappableType, Allocator>::erase_contained(const MappableType_& mappable)
{
    auto contained = this->contained_range(mappable);
    
    if (is_bidirectionally_sorted_ || ::size(contained) == bases(contained).size()) {
        this->erase(std::cbegin(contained).base(), std::cend(contained).base());
    } else {
        // TODO: find better implementation
        while (!contained.empty()) {
            this->erase(contained.front());
            contained = this->contained_range(mappable);
        }
    }
}

// non-member methods

template <typename MappableType, typename Allocator>
bool operator==(const MappableFlatSet<MappableType, Allocator>& lhs,
                const MappableFlatSet<MappableType, Allocator>& rhs)
{
    return lhs.elements_ == rhs.elements_;
}

template <typename MappableType, typename Allocator>
bool operator<(const MappableFlatSet<MappableType, Allocator>& lhs,
               const MappableFlatSet<MappableType, Allocator>& rhs)
{
    return lhs.elements_ < rhs.elements_;
}

template <typename MappableType, typename Allocator>
void swap(MappableFlatSet<MappableType, Allocator>& lhs,
          MappableFlatSet<MappableType, Allocator>& rhs)
{
    std::swap(lhs.elements_, rhs.elements_);
    std::swap(lhs.is_bidirectionally_sorted_, rhs.is_bidirectionally_sorted_);
    std::swap(lhs.max_element_size_, rhs.max_element_size_);
}

#endif /* mappable_flat_set_hpp */
