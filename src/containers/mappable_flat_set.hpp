// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

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

#include <concepts/comparable.hpp>
#include <concepts/mappable.hpp>
#include <concepts/mappable_range.hpp>
#include <utils/mappable_algorithms.hpp>

namespace octopus {

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
    
    MappableFlatSet(std::initializer_list<MappableType> mappables);
    
    MappableFlatSet(const MappableFlatSet&)            = default;
    MappableFlatSet& operator=(const MappableFlatSet&) = default;
    MappableFlatSet(MappableFlatSet&&)                 = default;
    MappableFlatSet& operator=(MappableFlatSet&&)      = default;
    
    ~MappableFlatSet() = default;
    
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
    std::pair<iterator, bool> emplace(Args...);
    std::pair<iterator, bool> insert(const MappableType&);
    std::pair<iterator, bool> insert(MappableType&&);
//    iterator insert(const_iterator, const MappableType& mappable);
//    iterator insert(const_iterator, MappableType&& mappable);
    template <typename InputIterator>
    void insert(InputIterator, InputIterator);
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
    
    iterator find(const MappableType&);
    const_iterator find(const MappableType&) const;
    size_type count(const MappableType&) const;
    
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
    typename RegionType<MappableType>::Position max_element_size_;
};

template <typename MappableType, typename Allocator>
MappableFlatSet<MappableType, Allocator>::MappableFlatSet()
:
elements_ {},
is_bidirectionally_sorted_ {true},
max_element_size_ {0}
{}

template <typename MappableType, typename Allocator>
template <typename InputIterator>
MappableFlatSet<MappableType, Allocator>::MappableFlatSet(InputIterator first, InputIterator second)
:
elements_ {first, second},
is_bidirectionally_sorted_ {true},
max_element_size_ {0}
{
    if (elements_.empty()) return;
    
    std::sort(std::begin(elements_), std::end(elements_));
    
    elements_.erase(std::unique(std::begin(elements_), std::end(elements_)), std::end(elements_));
    
    is_bidirectionally_sorted_ = is_bidirectionally_sorted(elements_);
    
    max_element_size_ = region_size(*largest_mappable(elements_));
}

template <typename MappableType, typename Allocator>
MappableFlatSet<MappableType, Allocator>::MappableFlatSet(std::initializer_list<MappableType> mappables)
:
elements_ {mappables},
is_bidirectionally_sorted_ {true},
max_element_size_ {0}
{
    if (elements_.empty()) return;
    
    std::sort(std::begin(elements_), std::end(elements_));
    
    elements_.erase(std::unique(std::begin(elements_), std::end(elements_)), std::end(elements_));
    
    is_bidirectionally_sorted_ = is_bidirectionally_sorted(elements_);
    
    max_element_size_ = region_size(*largest_mappable(elements_));
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

template <typename MappableType, typename Allocator>
template <typename ...Args>
std::pair<typename MappableFlatSet<MappableType, Allocator>::iterator, bool>
MappableFlatSet<MappableType, Allocator>::emplace(Args... args)
{
    elements_.emplace_back(std::forward<Args>(args)...);
    const auto it = std::lower_bound(std::begin(elements_), std::prev(std::end(elements_)),
                                     elements_.back());
    if (it != std::prev(std::end(elements_)) && *it == elements_.back()) {
        elements_.pop_back();
        return std::make_pair(it, false);
    }
    std::rotate(std::rbegin(elements_), std::next(std::rbegin(elements_)),
                std::make_reverse_iterator(it));
    if (is_bidirectionally_sorted_) {
        const auto overlapped = overlap_range(*it);
        is_bidirectionally_sorted_ = is_bidirectionally_sorted(overlapped);
    }
    max_element_size_ = std::max(max_element_size_, region_size(*it));
    return std::make_pair(it, true);
}

template <typename MappableType, typename Allocator>
std::pair<typename MappableFlatSet<MappableType, Allocator>::iterator, bool>
MappableFlatSet<MappableType, Allocator>::insert(const MappableType& m)
{
    auto it = std::lower_bound(std::begin(elements_), std::end(elements_), m);
    if (it == std::end(elements_) || !(*it == m)) {
        it = elements_.insert(it, m);
    } else {
        return std::make_pair(it, false);
    }
    if (is_bidirectionally_sorted_) {
        const auto overlapped = overlap_range(*it);
        is_bidirectionally_sorted_ = is_bidirectionally_sorted(overlapped);
    }
    max_element_size_ = std::max(max_element_size_, region_size(*it));
    return std::make_pair(it, true);
}

template <typename MappableType, typename Allocator>
std::pair<typename MappableFlatSet<MappableType, Allocator>::iterator, bool>
MappableFlatSet<MappableType, Allocator>::insert(MappableType&& m)
{
    auto it = std::lower_bound(std::begin(elements_), std::end(elements_), m);
    if (it == std::end(elements_) || !(*it == m)) {
        it = elements_.insert(it, m);
    } else {
        return std::make_pair(it, false);
    }
    if (is_bidirectionally_sorted_) {
        const auto overlapped = overlap_range(*it);
        is_bidirectionally_sorted_ = is_bidirectionally_sorted(overlapped);
    }
    max_element_size_ = std::max(max_element_size_, region_size(*it));
    return std::make_pair(it, true);
}

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

template <typename MappableType, typename Allocator>
template <typename InputIterator>
void MappableFlatSet<MappableType, Allocator>::insert(InputIterator first, InputIterator last)
{
    if (first == last) return;
    
    max_element_size_ = std::max(max_element_size_, region_size(*largest_mappable(first, last)));
    
    for (auto it1 = first; it1 != last; ) {
        const auto it2 = std::is_sorted_until(it1, last);
        
        auto ub = std::upper_bound(std::begin(elements_), std::end(elements_), *std::prev(it2));
        
        const auto d = std::distance(it1, it2);
        
        const auto it3 = elements_.insert(ub, it1, it2);
        
        // ub is now invalidated
        
        const auto lb = std::lower_bound(std::begin(elements_), it3, *it3);
        
        ub = std::next(it3, d);
        
        std::inplace_merge(lb, it3, ub);
        
        elements_.erase(std::unique(lb, ub), ub);
        
        it1 = it2;
    }
    
    if (is_bidirectionally_sorted_) {
        is_bidirectionally_sorted_ = is_bidirectionally_sorted(elements_);
    }
}

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
typename MappableFlatSet<MappableType, Allocator>::size_type
MappableFlatSet<MappableType, Allocator>::erase(const MappableType& m)
{
    const auto it = std::lower_bound(std::cbegin(elements_), std::cend(elements_), m);
    
    if (it != std::cend(elements_) && *it == m) {
        const auto m_size = region_size(m);
        
        elements_.erase(it);
        
        if (elements_.empty()) {
            max_element_size_ = 0;
            is_bidirectionally_sorted_ = true;
        } else {
            if (max_element_size_ == m_size) {
                max_element_size_ = region_size(*largest_mappable(elements_));
            }
        }
        
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
typename MappableFlatSet<MappableType, Allocator>::size_type
MappableFlatSet<MappableType, Allocator>::erase_all(InputIt first, const InputIt last)
{
    size_type num_erased {0};
    
    if (first == last) return num_erased;
    
    const auto region = encompassing_region(first, last);
    
    using octopus::contained_range;
    
    auto contained_elements = bases(contained_range(std::begin(elements_), std::end(elements_), region));
    
    if (contained_elements.empty()) return num_erased;
    
    typename RegionType<MappableType>::Size max_erased_size {0};
    
    auto first_contained = std::begin(contained_elements);
    auto last_contained  = std::end(contained_elements);
    
    auto last_element = std::end(elements_);
    
    while (first != last) {
        const auto it = std::lower_bound(first_contained, last_contained, *first);
        
        if (it != last_contained) {
            const auto p = std::mismatch(std::next(it), last_contained, std::next(first), last);
            
            const auto n = std::distance(p.first, last_contained);
            
            last_element = std::rotate(it, p.first, last_element);
            
            first_contained = it;
            last_contained  = std::next(it, n);
            
            const auto m = region_size(*largest_mappable(first, p.second));
            
            if (m > max_erased_size) {
                max_erased_size = m;
            }
            
            num_erased += std::distance(first, p.second);
            
            first = p.second;
        } else {
            ++first;
        }
    }
    
    if (num_erased > 0) {
        elements_.erase(last_element, std::end(elements_));
        
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
    
    return num_erased;
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
typename MappableFlatSet<MappableType, Allocator>::iterator
MappableFlatSet<MappableType, Allocator>::find(const MappableType& m)
{
    const auto it = std::lower_bound(std::begin(elements_), std::end(elements_), m);
    
    if (it == std::end(elements_) || !(*it == m)) return std::end(elements_);
    
    return it;
}

template <typename MappableType, typename Allocator>
typename MappableFlatSet<MappableType, Allocator>::const_iterator
MappableFlatSet<MappableType, Allocator>::find(const MappableType& m) const
{
    const auto it = std::lower_bound(std::cbegin(elements_), std::cend(elements_), m);
    
    if (it == std::cend(elements_) || !(*it == m)) return std::cend(elements_);
    
    return it;
}

template <typename MappableType, typename Allocator>
typename MappableFlatSet<MappableType, Allocator>::size_type
MappableFlatSet<MappableType, Allocator>::count(const MappableType& m) const
{
    return std::binary_search(std::cbegin(elements_), std::cend(elements_), m);
}

template <typename MappableType, typename Allocator>
const MappableType& MappableFlatSet<MappableType, Allocator>::leftmost() const
{
    return front();
}

template <typename MappableType, typename Allocator>
const MappableType& MappableFlatSet<MappableType, Allocator>::rightmost() const
{
    const auto& last = *std::prev(std::cend(elements_));
    
    if (is_bidirectionally_sorted_) {
        return last;
    } else {
        using octopus::overlap_range;
        const auto overlapped = overlap_range(std::cbegin(elements_), std::cend(elements_), last,
                                              max_element_size_);
        return *rightmost_mappable(std::cbegin(overlapped), std::cend(overlapped));
    }
}

template <typename MappableType, typename Allocator>
template <typename MappableType_>
bool
MappableFlatSet<MappableType, Allocator>::has_overlapped(const MappableType_& mappable) const
{
    using octopus::has_overlapped;
    if (is_bidirectionally_sorted_) {
        has_overlapped(std::cbegin(elements_), std::cend(elements_), mappable,
                       BidirectionallySortedTag {});
    }
    return has_overlapped(std::cbegin(elements_), std::cend(elements_), mappable);
}

template <typename MappableType, typename Allocator>
template <typename MappableType_>
bool
MappableFlatSet<MappableType, Allocator>::has_overlapped(iterator first, iterator last,
                                                         const MappableType_& mappable) const
{
    using octopus::has_overlapped;
    if (is_bidirectionally_sorted_) {
        has_overlapped(first, last, mappable, BidirectionallySortedTag {});
    }
    return has_overlapped(first, last, mappable, max_element_size_);
}

template <typename MappableType, typename Allocator>
template <typename MappableType_>
bool
MappableFlatSet<MappableType, Allocator>::has_overlapped(const_iterator first, const_iterator last,
                                                         const MappableType_& mappable) const
{
    using octopus::has_overlapped;
    if (is_bidirectionally_sorted_) {
        has_overlapped(first, last, mappable, BidirectionallySortedTag {});
    }
    return has_overlapped(first, last, mappable, max_element_size_);
}

template <typename MappableType, typename Allocator>
template <typename MappableType_>
typename MappableFlatSet<MappableType, Allocator>::size_type
MappableFlatSet<MappableType, Allocator>::count_overlapped(const MappableType_& mappable) const
{
    const auto overlapped = overlap_range(mappable);
    
    using octopus::size;
    
    if (is_bidirectionally_sorted_) {
        return size(overlapped, BidirectionallySortedTag {});
    }
    
    return size(overlapped);
}

template <typename MappableType, typename Allocator>
template <typename MappableType_>
typename MappableFlatSet<MappableType, Allocator>::size_type
MappableFlatSet<MappableType, Allocator>::count_overlapped(iterator first, iterator last,
                                                           const MappableType_& mappable) const
{
    const auto overlapped = overlap_range(first, last, mappable);
    
    using octopus::size;
    
    if (is_bidirectionally_sorted_) {
        return size(overlapped, BidirectionallySortedTag {});
    }
    
    return size(overlapped);
}

template <typename MappableType, typename Allocator>
template <typename MappableType_>
typename MappableFlatSet<MappableType, Allocator>::size_type
MappableFlatSet<MappableType, Allocator>::count_overlapped(const_iterator first, const_iterator last,
                                                           const MappableType_& mappable) const
{
    const auto overlapped = overlap_range(first, last, mappable);
    
    using octopus::size;
    
    if (is_bidirectionally_sorted_) {
        return size(overlapped, BidirectionallySortedTag {});
    }
    
    return size(overlapped);
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
    return overlap_range(const_iterator {first}, const_iterator {last}, mappable);
}

template <typename MappableType, typename Allocator>
template <typename MappableType_>
OverlapRange<typename MappableFlatSet<MappableType, Allocator>::const_iterator>
MappableFlatSet<MappableType, Allocator>::overlap_range(const_iterator first, const_iterator last,
                                                        const MappableType_& mappable) const
{
    using octopus::overlap_range;
    if (is_bidirectionally_sorted_) {
        overlap_range(first, last, mappable, BidirectionallySortedTag {});
    }
    return overlap_range(first, last, mappable, max_element_size_);
}

template <typename MappableType, typename Allocator>
template <typename MappableType_>
void MappableFlatSet<MappableType, Allocator>::erase_overlapped(const MappableType_& mappable)
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
    return has_contained(const_iterator {first}, const_iterator {last}, mappable);
}

template <typename MappableType, typename Allocator>
template <typename MappableType_>
bool
MappableFlatSet<MappableType, Allocator>::has_contained(const_iterator first, const_iterator last,
                                                             const MappableType_& mappable) const
{
    using octopus::has_contained;
    return has_contained(first, last, mappable);
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
    
    using octopus::size;
    
    if (is_bidirectionally_sorted_) {
        return size(contained, BidirectionallySortedTag {});
    }
    
    return size(contained);
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
    using octopus::contained_range;
    return contained_range(first, last, mappable);
}

template <typename MappableType, typename Allocator>
template <typename MappableType_>
void MappableFlatSet<MappableType, Allocator>::erase_contained(const MappableType_& mappable)
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

} // namespace octopus

#endif
