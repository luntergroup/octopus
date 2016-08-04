//
//  mappable_flat_multi_set.hpp
//  octopus
//
//  Created by Daniel Cooke on 19/07/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_mappable_set_hpp
#define Octopus_mappable_set_hpp

#include <memory>
#include <functional>
#include <algorithm>
#include <iterator>
#include <initializer_list>
#include <stdexcept>
#include <vector>

#include <boost/container/flat_set.hpp>

#include <concepts/comparable.hpp>
#include <concepts/mappable.hpp>
#include <concepts/mappable_range.hpp>
#include <utils/mappable_algorithms.hpp>

namespace octopus {

/*
 MappableFlatMultiSet is a container designed to allow fast retrival of MappableType elements with minimal
 memory overhead.
 */
template <typename MappableType, typename Allocator = std::allocator<MappableType>>
class MappableFlatMultiSet : public Comparable<MappableFlatMultiSet<MappableType, Allocator>>
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
    
    MappableFlatMultiSet();
    
    template <typename InputIterator>
    MappableFlatMultiSet(InputIterator first, InputIterator second);
    
    MappableFlatMultiSet(std::initializer_list<MappableType> mappables);
    
    MappableFlatMultiSet(const MappableFlatMultiSet&)            = default;
    MappableFlatMultiSet& operator=(const MappableFlatMultiSet&) = default;
    MappableFlatMultiSet(MappableFlatMultiSet&&)                 = default;
    MappableFlatMultiSet& operator=(MappableFlatMultiSet&&)      = default;
    
    ~MappableFlatMultiSet() = default;
    
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
    
    template <typename ...Args>
    iterator emplace(Args...);
    iterator insert(const MappableType&);
    iterator insert(MappableType&&);
    iterator insert(const_iterator, const MappableType& mappable);
    iterator insert(const_iterator, MappableType&& mappable);
    template <typename InputIterator>
    void insert(InputIterator, InputIterator);
    iterator insert(std::initializer_list<MappableType>);
    iterator erase(const_iterator);
    size_type erase(const MappableType&);
    iterator erase(const_iterator, const_iterator);
    template <typename InputIt>
    size_type erase_all(InputIt first, InputIt last);
    
    void clear();
    
    void swap(const MappableFlatMultiSet&);
    
    size_type size() const noexcept;
    size_type capacity() const noexcept;
    size_type max_size() const noexcept;
    bool empty() const noexcept;
    void reserve(size_type n);
    void shrink_to_fit();
    
    allocator_type get_allocator() noexcept;
    
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
    
    template <typename MappableType1_, typename MappableType2_>
    bool has_shared(const MappableType1_& mappable1, const MappableType2_& mappable2) const;
    template <typename MappableType1_, typename MappableType2_>
    bool has_shared(iterator first, iterator last,
                    const MappableType1_& mappable1, const MappableType2_& mappable2) const;
    template <typename MappableType1_, typename MappableType2_>
    bool has_shared(const_iterator first, const_iterator last,
                    const MappableType1_& mappable1, const MappableType2_& mappable2) const;
    
    template <typename MappableType1_, typename MappableType2_>
    size_type count_shared(const MappableType1_& mappable1, const MappableType2_& mappable2) const;
    template <typename MappableType1_, typename MappableType2_>
    size_type count_shared(iterator first, iterator last,
                           const MappableType1_& mappable1, const MappableType2_& mappable2) const;
    template <typename MappableType1_, typename MappableType2_>
    size_type count_shared(const_iterator first, const_iterator last,
                           const MappableType1_& mappable1, const MappableType2_& mappable2) const;
    
    template <typename MappableType1_, typename MappableType2_>
    SharedRange<const_iterator> shared_range(const MappableType1_& mappable1,
                                             const MappableType2_& mappable2) const;
    template <typename MappableType1_, typename MappableType2_>
    SharedRange<iterator> shared_range(iterator first, iterator last,
                                       const MappableType1_& mappable1,
                                       const MappableType2_& mappable2) const;
    template <typename MappableType1_, typename MappableType2_>
    SharedRange<const_iterator> shared_range(const_iterator first, const_iterator last,
                                             const MappableType1_& mappable1,
                                             const MappableType2_& mappable2) const;
    
    template <typename M, typename A>
    friend bool operator==(const MappableFlatMultiSet<M, A>& lhs, const MappableFlatMultiSet<M, A>& rhs);
    template <typename M, typename A>
    friend bool operator<(const MappableFlatMultiSet<M, A>& lhs, const MappableFlatMultiSet<M, A>& rhs);
    template <typename M, typename A>
    friend void swap(MappableFlatMultiSet<M, A>& lhs, MappableFlatMultiSet<M, A>& rhs);
    
private:
    base_t elements_;
    bool is_bidirectionally_sorted_;
    typename RegionType<MappableType>::Position max_element_size_;
};

template <typename MappableType, typename Allocator>
MappableFlatMultiSet<MappableType, Allocator>::MappableFlatMultiSet()
:
elements_ {},
is_bidirectionally_sorted_ {true},
max_element_size_ {}
{}

template <typename MappableType, typename Allocator>
template <typename InputIterator>
MappableFlatMultiSet<MappableType, Allocator>::MappableFlatMultiSet(InputIterator first, InputIterator second)
:
elements_ {first, second},
is_bidirectionally_sorted_ {is_bidirectionally_sorted(elements_)},
max_element_size_ {(elements_.empty()) ? 0 : region_size(*largest_mappable(elements_))}
{}

template <typename MappableType, typename Allocator>
MappableFlatMultiSet<MappableType, Allocator>::MappableFlatMultiSet(std::initializer_list<MappableType> mappables)
:
elements_ {mappables},
is_bidirectionally_sorted_ {is_bidirectionally_sorted(elements_)},
max_element_size_ {(elements_.empty()) ? 0 : region_size(*largest_mappable(elements_))}
{}

template <typename MappableType, typename Allocator>
typename MappableFlatMultiSet<MappableType, Allocator>::iterator
MappableFlatMultiSet<MappableType, Allocator>::begin() noexcept
{
    return elements_.begin();
}

template <typename MappableType, typename Allocator>
typename MappableFlatMultiSet<MappableType, Allocator>::const_iterator
MappableFlatMultiSet<MappableType, Allocator>::begin() const noexcept
{
    return elements_.begin();
}

template <typename MappableType, typename Allocator>
typename MappableFlatMultiSet<MappableType, Allocator>::const_iterator
MappableFlatMultiSet<MappableType, Allocator>::cbegin() const noexcept
{
    return elements_.cbegin();
}

template <typename MappableType, typename Allocator>
typename MappableFlatMultiSet<MappableType, Allocator>::iterator
MappableFlatMultiSet<MappableType, Allocator>::end() noexcept
{
    return elements_.end();
}

template <typename MappableType, typename Allocator>
typename MappableFlatMultiSet<MappableType, Allocator>::const_iterator
MappableFlatMultiSet<MappableType, Allocator>::end() const noexcept
{
    return elements_.end();
}

template <typename MappableType, typename Allocator>
typename MappableFlatMultiSet<MappableType, Allocator>::const_iterator
MappableFlatMultiSet<MappableType, Allocator>::cend() const noexcept
{
    return elements_.cend();
}

template <typename MappableType, typename Allocator>
typename MappableFlatMultiSet<MappableType, Allocator>::reverse_iterator
MappableFlatMultiSet<MappableType, Allocator>::rbegin() noexcept
{
    return elements_.rbegin();
}

template <typename MappableType, typename Allocator>
typename MappableFlatMultiSet<MappableType, Allocator>::const_reverse_iterator
MappableFlatMultiSet<MappableType, Allocator>::rbegin() const noexcept
{
    return elements_.rbegin();
}

template <typename MappableType, typename Allocator>
typename MappableFlatMultiSet<MappableType, Allocator>::const_reverse_iterator
MappableFlatMultiSet<MappableType, Allocator>::crbegin() const noexcept
{
    return elements_.crbegin();
}

template <typename MappableType, typename Allocator>
typename MappableFlatMultiSet<MappableType, Allocator>::reverse_iterator
MappableFlatMultiSet<MappableType, Allocator>::rend() noexcept
{
    return elements_.rend();
}

template <typename MappableType, typename Allocator>
typename MappableFlatMultiSet<MappableType, Allocator>::const_reverse_iterator
MappableFlatMultiSet<MappableType, Allocator>::rend() const noexcept
{
    return elements_.rend();
}

template <typename MappableType, typename Allocator>
typename MappableFlatMultiSet<MappableType, Allocator>::const_reverse_iterator
MappableFlatMultiSet<MappableType, Allocator>::crend() const noexcept
{
    return elements_.crend();
}

template <typename MappableType, typename Allocator>
typename MappableFlatMultiSet<MappableType, Allocator>::reference
MappableFlatMultiSet<MappableType, Allocator>::at(size_type pos)
{
    if (pos < size()) {
        return *std::next(begin(), pos);
    } else {
        throw std::out_of_range {"MappableFlatMultiSet"};
    }
}

template <typename MappableType, typename Allocator>
typename MappableFlatMultiSet<MappableType, Allocator>::const_reference
MappableFlatMultiSet<MappableType, Allocator>::at(size_type pos) const
{
    if (pos < size()) {
        return *std::next(cbegin(), pos);
    } else {
        throw std::out_of_range {"MappableFlatMultiSet"};
    }
}

template <typename MappableType, typename Allocator>
typename MappableFlatMultiSet<MappableType, Allocator>::reference
MappableFlatMultiSet<MappableType, Allocator>::operator[](size_type pos)
{
    return *std::next(begin(), pos);
}

template <typename MappableType, typename Allocator>
typename MappableFlatMultiSet<MappableType, Allocator>::const_reference
MappableFlatMultiSet<MappableType, Allocator>::operator[](size_type pos) const
{
    return *std::next(cbegin(), pos);
}

template <typename MappableType, typename Allocator>
typename MappableFlatMultiSet<MappableType, Allocator>::reference
MappableFlatMultiSet<MappableType, Allocator>::front()
{
    return *begin();
}

template <typename MappableType, typename Allocator>
typename MappableFlatMultiSet<MappableType, Allocator>::const_reference
MappableFlatMultiSet<MappableType, Allocator>::front() const
{
    return *cbegin();
}

template <typename MappableType, typename Allocator>
typename MappableFlatMultiSet<MappableType, Allocator>::reference
MappableFlatMultiSet<MappableType, Allocator>::back()
{
    return *std::prev(end());
}

template <typename MappableType, typename Allocator>
typename MappableFlatMultiSet<MappableType, Allocator>::const_reference
MappableFlatMultiSet<MappableType, Allocator>::back() const
{
    return *std::prev(cend());
}

template <typename MappableType, typename Allocator>
template <typename ...Args>
typename MappableFlatMultiSet<MappableType, Allocator>::iterator
MappableFlatMultiSet<MappableType, Allocator>::emplace(Args... args)
{
    const auto it = elements_.emplace(std::forward<Args>(args)...);
    if (is_bidirectionally_sorted_) {
        const auto overlapped = overlap_range(*it);
        is_bidirectionally_sorted_ = is_bidirectionally_sorted(overlapped);
    }
    max_element_size_ = std::max(max_element_size_, region_size(*it));
    return it;
}

template <typename MappableType, typename Allocator>
typename MappableFlatMultiSet<MappableType, Allocator>::iterator
MappableFlatMultiSet<MappableType, Allocator>::insert(const MappableType& m)
{
    const auto it = elements_.insert(m);
    if (is_bidirectionally_sorted_) {
        const auto overlapped = overlap_range(*it);
        is_bidirectionally_sorted_ = is_bidirectionally_sorted(overlapped);
    }
    max_element_size_ = std::max(max_element_size_, region_size(*it));
    return it;
}

template <typename MappableType, typename Allocator>
typename MappableFlatMultiSet<MappableType, Allocator>::iterator
MappableFlatMultiSet<MappableType, Allocator>::insert(MappableType&& m)
{
    const auto it = elements_.insert(std::move(m));
    if (is_bidirectionally_sorted_) {
        const auto overlapped = overlap_range(*it);
        is_bidirectionally_sorted_ = is_bidirectionally_sorted(overlapped);
    }
    max_element_size_ = std::max(max_element_size_, region_size(*it));
    return it;
}

template <typename MappableType, typename Allocator>
typename MappableFlatMultiSet<MappableType, Allocator>::iterator
MappableFlatMultiSet<MappableType, Allocator>::insert(const_iterator it, const MappableType& m)
{
    const auto it2 = elements_.insert(it, m);
    if (is_bidirectionally_sorted_) {
        const auto overlapped = overlap_range(*it2);
        is_bidirectionally_sorted_ = is_bidirectionally_sorted(overlapped);
    }
    max_element_size_ = std::max(max_element_size_, region_size(*it));
    return it2;
}

template <typename MappableType, typename Allocator>
typename MappableFlatMultiSet<MappableType, Allocator>::iterator
MappableFlatMultiSet<MappableType, Allocator>::insert(const_iterator it, MappableType&& m)
{
    const auto it2 = elements_.insert(it, std::move(m));
    if (is_bidirectionally_sorted_) {
        const auto overlapped = overlap_range(*it2);
        is_bidirectionally_sorted_ = is_bidirectionally_sorted(overlapped);
    }
    max_element_size_ = std::max(max_element_size_, region_size(*it));
    return it2;
}

template <typename MappableType, typename Allocator>
template <typename InputIterator>
void
MappableFlatMultiSet<MappableType, Allocator>::insert(InputIterator first, InputIterator last)
{
    if (first != last) {
        max_element_size_ = std::max(max_element_size_, region_size(*largest_mappable(first, last)));
    }
    elements_.insert(first, last);
    if (is_bidirectionally_sorted_) {
        is_bidirectionally_sorted_ = is_bidirectionally_sorted(elements_);
    }
}

template <typename MappableType, typename Allocator>
typename MappableFlatMultiSet<MappableType, Allocator>::iterator
MappableFlatMultiSet<MappableType, Allocator>::insert(std::initializer_list<MappableType> il)
{
    if (!il.empty()) {
        max_element_size_ = std::max(max_element_size_, region_size(*largest_element(il)));
    }
    return elements_.insert(std::move(il));
    if (is_bidirectionally_sorted_) {
        is_bidirectionally_sorted_ = is_bidirectionally_sorted(elements_);
    }
}

template <typename MappableType, typename Allocator>
typename MappableFlatMultiSet<MappableType, Allocator>::iterator
MappableFlatMultiSet<MappableType, Allocator>::erase(const_iterator p)
{
    if (p == cend()) return elements_.erase(p);
    
    const auto erased_size = region_size(*p);
    
    const auto result = elements_.erase(p);
    
    if (elements_.empty()) {
        max_element_size_ = 0;
        is_bidirectionally_sorted_ = true;
    } else {
        if (max_element_size_ == erased_size) {
            max_element_size_ = region_size(*largest_mappable(elements_));
        }
    }
    
    return result;
}

template <typename MappableType, typename Allocator>
typename MappableFlatMultiSet<MappableType, Allocator>::size_type
MappableFlatMultiSet<MappableType, Allocator>::erase(const MappableType& m)
{
    const auto m_size = region_size(m);
    
    const auto result = elements_.erase(m);
    
    if (result > 0) {
        if (elements_.empty()) {
            max_element_size_ = 0;
            is_bidirectionally_sorted_ = true;
        } else {
            if (max_element_size_ == m_size) {
                max_element_size_ = region_size(*largest_mappable(elements_));
            }
        }
        
        return result;
    }
    
    return 0;
}

template <typename MappableType, typename Allocator>
typename MappableFlatMultiSet<MappableType, Allocator>::iterator
MappableFlatMultiSet<MappableType, Allocator>::erase(const_iterator first, const_iterator last)
{
    if (first == last) return elements_.erase(first, last);
    
    const auto max_erased_size = region_size(*largest_mappable(first, last));
    
    const auto result = elements_.erase(first, last);
    
    if (elements_.empty()) {
        max_element_size_ = 0;
        is_bidirectionally_sorted_ = true;
    } else {
        if (max_element_size_ == max_erased_size) {
            max_element_size_ = region_size(*largest_mappable(elements_));
        }
    }
    
    return result;
}

template <typename MappableType, typename Allocator>
template <typename InputIt>
typename MappableFlatMultiSet<MappableType, Allocator>::size_type
MappableFlatMultiSet<MappableType, Allocator>::erase_all(InputIt first, InputIt last)
{
    size_type result {0};
    
    if (first == last) return result;
    
    auto from = std::cbegin(elements_);
    
    typename RegionType<MappableType>::SizeType max_erased_size {0};
    
    std::for_each(first, last, [this, &result, &from, &max_erased_size] (const auto& element) {
        const auto er = elements_.equal_range(element);
        
        if (er.first != er.second) {
            if (region_size(element) > max_erased_size) {
                max_erased_size = region_size(element);
            }
            
            result += std::distance(er.first, er.second);
            
            elements_.erase(er.first, er.second);
        }
    });
    
    if (result > 0) {
        if (!elements_.empty()) {
            if (!is_bidirectionally_sorted_) {
                is_bidirectionally_sorted_ = is_bidirectionally_sorted(elements_);
            }
            if (max_element_size_ == max_erased_size) {
                max_element_size_ = region_size(*largest_mappable(elements_));
            }
        } else {
            max_element_size_ = 0;
            is_bidirectionally_sorted_ = true;
        }
    }
    
    return result;
}

template <typename MappableType, typename Allocator>
void MappableFlatMultiSet<MappableType, Allocator>::clear()
{
    elements_.clear();
    is_bidirectionally_sorted_ = true;
    max_element_size_ = 0;
}

template <typename MappableType, typename Allocator>
void MappableFlatMultiSet<MappableType, Allocator>::swap(const MappableFlatMultiSet& m)
{
    std::swap(elements_, m.elements_);
    std::swap(is_bidirectionally_sorted_, m.is_bidirectionally_sorted_);
    std::swap(max_element_size_, m.max_element_size_);
}

template <typename MappableType, typename Allocator>
typename MappableFlatMultiSet<MappableType, Allocator>::size_type
MappableFlatMultiSet<MappableType, Allocator>::size() const noexcept
{
    return elements_.size();
}

template <typename MappableType, typename Allocator>
typename MappableFlatMultiSet<MappableType, Allocator>::size_type
MappableFlatMultiSet<MappableType, Allocator>::capacity() const noexcept
{
    return elements_.capacity();
}

template <typename MappableType, typename Allocator>
typename MappableFlatMultiSet<MappableType, Allocator>::size_type
MappableFlatMultiSet<MappableType, Allocator>::max_size() const noexcept
{
    return elements_.max_size();
}

template <typename MappableType, typename Allocator>
bool MappableFlatMultiSet<MappableType, Allocator>::empty() const noexcept
{
    return elements_.empty();
}

template <typename MappableType, typename Allocator>
void
MappableFlatMultiSet<MappableType, Allocator>::reserve(size_type n)
{
    elements_.reserve(n);
}

template <typename MappableType, typename Allocator>
void
MappableFlatMultiSet<MappableType, Allocator>::shrink_to_fit()
{
    elements_.shrink_to_fit();
}

template <typename MappableType, typename Allocator>
typename MappableFlatMultiSet<MappableType, Allocator>::allocator_type
MappableFlatMultiSet<MappableType, Allocator>::get_allocator() noexcept
{
    return elements_.get_allocator();
}

template <typename MappableType, typename Allocator>
const MappableType& MappableFlatMultiSet<MappableType, Allocator>::leftmost() const
{
    return front();
}

template <typename MappableType, typename Allocator>
const MappableType& MappableFlatMultiSet<MappableType, Allocator>::rightmost() const
{
    using std::cbegin; using std::cend;
    const auto& last = *std::prev(elements_.cend());
    if (is_bidirectionally_sorted_) {
        return last;
    } else {
        using octopus::overlap_range;
        const auto overlapped = overlap_range(cbegin(elements_), cend(elements_), last,
                                              max_element_size_);
        return *rightmost_mappable(cbegin(overlapped), cend(overlapped));
    }
}

template <typename MappableType, typename Allocator>
template <typename MappableType_>
bool
MappableFlatMultiSet<MappableType, Allocator>::has_overlapped(const MappableType_& mappable) const
{
    using octopus::has_overlapped;
    if (is_bidirectionally_sorted_) {
        has_overlapped(std::begin(elements_), std::end(elements_), mappable, BidirectionallySortedTag {});
    }
    return has_overlapped(std::begin(elements_), std::end(elements_), mappable);
}

template <typename MappableType, typename Allocator>
template <typename MappableType_>
bool
MappableFlatMultiSet<MappableType, Allocator>::has_overlapped(iterator first, iterator last,
                                                              const MappableType_& mappable) const
{
    using octopus::has_overlapped;
    if (is_bidirectionally_sorted_) {
        has_overlapped(first, last, mappable, BidirectionallySortedTag {});
    }
    return has_overlapped(first, last, mappable);
}

template <typename MappableType, typename Allocator>
template <typename MappableType_>
bool
MappableFlatMultiSet<MappableType, Allocator>::has_overlapped(const_iterator first, const_iterator last,
                                                              const MappableType_& mappable) const
{
    using octopus::has_overlapped;
    if (is_bidirectionally_sorted_) {
        has_overlapped(first, last, mappable, BidirectionallySortedTag {});
    }
    return has_overlapped(first, last, mappable);
}

template <typename MappableType, typename Allocator>
template <typename MappableType_>
typename MappableFlatMultiSet<MappableType, Allocator>::size_type
MappableFlatMultiSet<MappableType, Allocator>::count_overlapped(const MappableType_& mappable) const
{
    using octopus::size;
    
    const auto overlapped = overlap_range(mappable);
    
    if (is_bidirectionally_sorted_) {
        return size(overlapped, BidirectionallySortedTag {});
    }
    
    return size(overlapped);
}

template <typename MappableType, typename Allocator>
template <typename MappableType_>
typename MappableFlatMultiSet<MappableType, Allocator>::size_type
MappableFlatMultiSet<MappableType, Allocator>::count_overlapped(iterator first, iterator last,
                                                                const MappableType_& mappable) const
{
    using octopus::size;
    
    const auto overlapped = overlap_range(first, last, mappable);
    
    if (is_bidirectionally_sorted_) {
        return size(overlapped, BidirectionallySortedTag {});
    }
    
    return size(overlapped);
}

template <typename MappableType, typename Allocator>
template <typename MappableType_>
typename MappableFlatMultiSet<MappableType, Allocator>::size_type
MappableFlatMultiSet<MappableType, Allocator>::count_overlapped(const_iterator first, const_iterator last,
                                                                const MappableType_& mappable) const
{
    using octopus::size;
    
    const auto overlapped = overlap_range(first, last, mappable);
    
    if (is_bidirectionally_sorted_) {
        return size(overlapped, BidirectionallySortedTag {});
    }
    
    return size(overlapped);
}

template <typename MappableType, typename Allocator>
template <typename MappableType_>
OverlapRange<typename MappableFlatMultiSet<MappableType, Allocator>::const_iterator>
MappableFlatMultiSet<MappableType, Allocator>::overlap_range(const MappableType_& mappable) const
{
    return overlap_range(std::cbegin(elements_), std::cend(elements_), mappable);
}

template <typename MappableType, typename Allocator>
template <typename MappableType_>
OverlapRange<typename MappableFlatMultiSet<MappableType, Allocator>::iterator>
MappableFlatMultiSet<MappableType, Allocator>::overlap_range(iterator first, iterator last,
                                                             const MappableType_& mappable) const
{
    return overlap_range(const_iterator(first), const_iterator(last), mappable);
}

template <typename MappableType, typename Allocator>
template <typename MappableType_>
OverlapRange<typename MappableFlatMultiSet<MappableType, Allocator>::const_iterator>
MappableFlatMultiSet<MappableType, Allocator>::overlap_range(const_iterator first, const_iterator last,
                                                             const MappableType_& mappable) const
{
    using octopus::overlap_range;
    if (is_bidirectionally_sorted_) {
        return overlap_range(first, last, mappable, BidirectionallySortedTag {});
    }
    return overlap_range(first, last, mappable, max_element_size_);
}

template <typename MappableType, typename Allocator>
template <typename MappableType_>
void MappableFlatMultiSet<MappableType, Allocator>::erase_overlapped(const MappableType_& mappable)
{
    // TODO: find better implementation
    
    const auto overlapped = overlap_range(mappable);
    
    using octopus::size;
    
    if (is_bidirectionally_sorted_ || size(overlapped) == bases(overlapped).size()) {
        erase(std::cbegin(overlapped).base(), std::cend(overlapped).base());
    } else {
        const std::vector<MappableType> tmp {std::cbegin(overlapped), std::cend(overlapped)};
        
        erase_all(std::cbegin(tmp), std::cend(tmp));
    }
}

template <typename MappableType, typename Allocator>
template <typename MappableType_>
bool
MappableFlatMultiSet<MappableType, Allocator>::has_contained(const MappableType_& mappable) const
{
    return has_contained(std::cbegin(elements_), std::cend(elements_), mappable);
}

template <typename MappableType, typename Allocator>
template <typename MappableType_>
bool
MappableFlatMultiSet<MappableType, Allocator>::has_contained(iterator first, iterator last,
                                                             const MappableType_& mappable) const
{
    return has_contained(const_iterator(first), const_iterator(last), mappable);
}

template <typename MappableType, typename Allocator>
template <typename MappableType_>
bool
MappableFlatMultiSet<MappableType, Allocator>::has_contained(const_iterator first, const_iterator last,
                                                             const MappableType_& mappable) const
{
    using octopus::has_contained;
    return has_contained(first, last, mappable);
}

template <typename MappableType, typename Allocator>
template <typename MappableType_>
typename MappableFlatMultiSet<MappableType, Allocator>::size_type
MappableFlatMultiSet<MappableType, Allocator>::count_contained(const MappableType_& mappable) const
{
    return count_contained(std::cbegin(elements_), std::cend(elements_), mappable);
}

template <typename MappableType, typename Allocator>
template <typename MappableType_>
typename MappableFlatMultiSet<MappableType, Allocator>::size_type
MappableFlatMultiSet<MappableType, Allocator>::count_contained(iterator first, iterator last,
                                                               const MappableType_& mappable) const
{
    return count_contained(const_iterator(first), const_iterator(last), mappable);
}

template <typename MappableType, typename Allocator>
template <typename MappableType_>
typename MappableFlatMultiSet<MappableType, Allocator>::size_type
MappableFlatMultiSet<MappableType, Allocator>::count_contained(const_iterator first, const_iterator last,
                                                               const MappableType_& mappable) const
{
    const auto contained = contained_range(first, last, mappable);
    
    using octopus::size;
    
    if (is_bidirectionally_sorted_) {
        return size(contained, BidirectionallySortedTag {});
    }
    
    return size(contained);
}

template <typename MappableType, typename Allocator>
template <typename MappableType_>
ContainedRange<typename MappableFlatMultiSet<MappableType, Allocator>::const_iterator>
MappableFlatMultiSet<MappableType, Allocator>::contained_range(const MappableType_& mappable) const
{
    return contained_range(std::cbegin(elements_), std::cend(elements_), mappable);
}

template <typename MappableType, typename Allocator>
template <typename MappableType_>
ContainedRange<typename MappableFlatMultiSet<MappableType, Allocator>::iterator>
MappableFlatMultiSet<MappableType, Allocator>::contained_range(iterator first, iterator last,
                                                               const MappableType_& mappable) const
{
    return contained_range(const_iterator(first), const_iterator(last), mappable);
}

template <typename MappableType, typename Allocator>
template <typename MappableType_>
ContainedRange<typename MappableFlatMultiSet<MappableType, Allocator>::const_iterator>
MappableFlatMultiSet<MappableType, Allocator>::contained_range(const_iterator first, const_iterator last,
                                                               const MappableType_& mappable) const
{
    using octopus::contained_range;
    return contained_range(first, last, mappable);
}

template <typename MappableType, typename Allocator>
template <typename MappableType_>
void MappableFlatMultiSet<MappableType, Allocator>::erase_contained(const MappableType_& mappable)
{
    auto contained = this->contained_range(mappable);
    
    using octopus::size;
    
    if (is_bidirectionally_sorted_ || size(contained) == bases(contained).size()) {
        this->erase(std::cbegin(contained).base(), std::cend(contained).base());
    } else {
        // TODO: find better implementation
        while (!contained.empty()) {
            this->erase(contained.front());
            contained = this->contained_range(mappable);
        }
    }
}

template <typename MappableType, typename Allocator>
template <typename MappableType1_, typename MappableType2_>
bool
MappableFlatMultiSet<MappableType, Allocator>::has_shared(const MappableType1_& mappable1,
                                                          const MappableType2_& mappable2) const
{
    return has_shared(std::cbegin(elements_), std::cend(elements_), mappable1, mappable2);
}

template <typename MappableType, typename Allocator>
template <typename MappableType1_, typename MappableType2_>
bool
MappableFlatMultiSet<MappableType, Allocator>::has_shared(iterator first, iterator last,
                                                          const MappableType1_& mappable1,
                                                          const MappableType2_& mappable2) const
{
    return has_shared(const_iterator(first), const_iterator(last), mappable1, mappable2);
}

template <typename MappableType, typename Allocator>
template <typename MappableType1_, typename MappableType2_>
bool
MappableFlatMultiSet<MappableType, Allocator>::has_shared(const_iterator first, const_iterator last,
                                                          const MappableType1_& mappable1,
                                                          const MappableType2_& mappable2) const
{
    if (inner_distance(mappable1, mappable2) > max_element_size_) return false;
    
    const auto m = std::minmax(mapped_region(mappable1), mapped_region(mappable2));
    
    const auto overlapped_lhs = overlap_range(first, last, m.first);
    
    return std::any_of(std::cbegin(overlapped_lhs), std::cend(overlapped_lhs),
                       [&m] (const auto& region) { return overlaps(region, m.second); });
}

template <typename MappableType, typename Allocator>
template <typename MappableType1_, typename MappableType2_>
typename MappableFlatMultiSet<MappableType, Allocator>::size_type
MappableFlatMultiSet<MappableType, Allocator>::count_shared(const MappableType1_& mappable1,
                                                            const MappableType2_& mappable2) const
{
    return count_shared(std::cbegin(elements_), std::cend(elements_), mappable1, mappable2);
}

template <typename MappableType, typename Allocator>
template <typename MappableType1_, typename MappableType2_>
typename MappableFlatMultiSet<MappableType, Allocator>::size_type
MappableFlatMultiSet<MappableType, Allocator>::count_shared(iterator first, iterator last,
                                                            const MappableType1_& mappable1,
                                                            const MappableType2_& mappable2) const
{
    return count_shared(const_iterator(first), const_iterator(last), mappable1, mappable2);
}

template <typename MappableType, typename Allocator>
template <typename MappableType1_, typename MappableType2_>
typename MappableFlatMultiSet<MappableType, Allocator>::size_type
MappableFlatMultiSet<MappableType, Allocator>::count_shared(const_iterator first, const_iterator last,
                                                            const MappableType1_& mappable1,
                                                            const MappableType2_& mappable2) const
{
    if (inner_distance(mappable1, mappable2) > max_element_size_) return 0;
    
    const auto m = std::minmax(mapped_region(mappable1), mapped_region(mappable2));
    
    const auto overlapped_lhs = overlap_range(first, last, m.first);
    
    return std::count_if(std::cbegin(overlapped_lhs), std::cend(overlapped_lhs),
                       [&m] (const auto& region) { return overlaps(region, m.second); });
}

template <typename MappableType, typename Allocator>
template <typename MappableType1_, typename MappableType2_>
SharedRange<typename MappableFlatMultiSet<MappableType, Allocator>::const_iterator>
MappableFlatMultiSet<MappableType, Allocator>::shared_range(const MappableType1_& mappable1,
                                                            const MappableType2_& mappable2) const
{
    return shared_range(std::cbegin(elements_), std::cend(elements_), mappable1, mappable2);
}

template <typename MappableType, typename Allocator>
template <typename MappableType1_, typename MappableType2_>
SharedRange<typename MappableFlatMultiSet<MappableType, Allocator>::iterator>
MappableFlatMultiSet<MappableType, Allocator>::shared_range(iterator first, iterator last,
                                                            const MappableType1_& mappable1,
                                                            const MappableType2_& mappable2) const
{
    return shared_range(const_iterator(first), const_iterator(last), mappable1, mappable2);
}

template <typename MappableType, typename Allocator>
template <typename MappableType1_, typename MappableType2_>
SharedRange<typename MappableFlatMultiSet<MappableType, Allocator>::const_iterator>
MappableFlatMultiSet<MappableType, Allocator>::shared_range(const_iterator first, const_iterator last,
                                                            const MappableType1_& mappable1,
                                                            const MappableType2_& mappable2) const
{
    if (inner_distance(mappable1, mappable2) > max_element_size_) {
        return make_shared_range(last, last, mappable1, mappable2);
    }
    
    const auto m = std::minmax(mapped_region(mappable1), mapped_region(mappable2));
    
    const auto overlapped_lhs = overlap_range(first, last, m.first);
    
    const auto it = std::find_if(std::cbegin(overlapped_lhs), std::cend(overlapped_lhs),
                                 [&m] (const auto& region) { return overlaps(region, m.second); });
    
    auto end = std::prev(overlapped_lhs.end());
    
    while (end != it && !overlaps(*end, m.second)) --end;
    
    return make_shared_range(it.base(), std::next(end).base(), mappable1, mappable2);
}

// non-member methods

template <typename MappableType, typename Allocator>
bool operator==(const MappableFlatMultiSet<MappableType, Allocator>& lhs,
                const MappableFlatMultiSet<MappableType, Allocator>& rhs)
{
    return lhs.elements_ == rhs.elements_;
}

template <typename MappableType, typename Allocator>
bool operator<(const MappableFlatMultiSet<MappableType, Allocator>& lhs,
               const MappableFlatMultiSet<MappableType, Allocator>& rhs)
{
    return lhs.elements_ < rhs.elements_;
}

template <typename MappableType, typename Allocator>
void swap(MappableFlatMultiSet<MappableType, Allocator>& lhs,
          MappableFlatMultiSet<MappableType, Allocator>& rhs)
{
    std::swap(lhs.elements_, rhs.elements_);
    std::swap(lhs.is_bidirectionally_sorted_, rhs.is_bidirectionally_sorted_);
    std::swap(lhs.max_element_size_, rhs.max_element_size_);
}

template <typename ForwardIterator, typename MappableType1, typename MappableType2, typename Allocator>
ForwardIterator
find_first_shared(const MappableFlatMultiSet<MappableType1, Allocator>& mappables,
                  ForwardIterator first, ForwardIterator last,
                  const MappableType2& mappable)
{
    return std::find_if(first, last, [&mappables, &mappable] (const auto& m) {
                            return mappables.has_shared(m, mappable);
                        });
}

template <typename MappableType>
auto calculate_positional_coverage(const MappableFlatMultiSet<MappableType>& mappables,
                                   const GenomicRegion& region)
{
    const auto overlapped = mappables.overlap_range(region);
    return calculate_positional_coverage(overlapped, region);
}

} // namespace octopus

#endif
