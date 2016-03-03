//
//  mappable_map.hpp
//  Octopus
//
//  Created by Daniel Cooke on 20/07/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_mappable_map_hpp
#define Octopus_mappable_map_hpp

#include <unordered_map>
#include <algorithm>
#include <numeric>
#include <iterator>
#include <stdexcept>
#include <cstddef>
#include <functional>
#include <utility>

#include "mappable_flat_multi_set.hpp"

template <typename KeyType, typename MappableType, typename Allocator = std::allocator<MappableType>>
using MappableMap = std::unordered_map<KeyType, MappableFlatMultiSet<MappableType, Allocator>>;

template <typename Map>
MappableMap<typename Map::key_type, typename Map::mapped_type::value_type>
make_mappable_map(Map map)
{
    using std::make_move_iterator; using std::begin; using std::end;
    
    using MappableTp = typename Map::mapped_type::value_type;
    
    MappableMap<typename Map::key_type, MappableTp> result {map.size()};
    
    std::for_each(make_move_iterator(begin(map)), make_move_iterator(end(map)),
                  [&result] (auto&& p) {
                      result.emplace(std::piecewise_construct,
                                     std::forward_as_tuple(std::move(p.first)),
                                     std::forward_as_tuple(make_move_iterator(begin(p.second)),
                                                           make_move_iterator(end(p.second))));
                  });
    
    return result;
}

template <typename KeyType, typename MappableType>
auto encompassing_region(const MappableMap<KeyType, MappableType>& mappables)
{
    if (mappables.empty()) {
        throw std::runtime_error {"get_encompassing_region called with empty MappableMap"};
    }
    
    auto result = encompassing_region(std::cbegin(mappables)->second);
    
    std::for_each(std::next(std::cbegin(mappables)), std::cend(mappables),
                  [&result] (const auto& p) {
                      result = encompassing_region(result, encompassing_region(p.second));
                  });
    
    return result;
}

template <typename KeyType, typename MappableType>
size_t count_mappables(const MappableMap<KeyType, MappableType>& mappables)
{
    return std::accumulate(std::cbegin(mappables), std::cend(mappables), size_t {},
                           [] (const auto prev, const auto& v) { return prev + v.second.size(); });
}

template <typename KeyType, typename MappableType1, typename MappableType2>
bool
has_overlapped(const MappableMap<KeyType, MappableType1>& mappables, const MappableType2& mappable)
{
    return std::any_of(std::cbegin(mappables), std::cend(mappables),
                       [&mappable] (const auto& p) {
                           return p.second.has_overlapped(mappable);
                       });
}

template <typename KeyType, typename MappableType1, typename MappableType2>
typename MappableFlatMultiSet<MappableType1>::size_type
count_overlapped(const MappableMap<KeyType, MappableType1>& mappables, const MappableType2& mappable)
{
    using SizeType = typename MappableFlatMultiSet<MappableType1>::size_type;
    
    return std::accumulate(std::cbegin(mappables), std::cend(mappables), SizeType {0},
                           [&mappable] (const auto curr, const auto& p) {
                               return curr + p.second.count_overlapped(mappable);
                           });
}

template <typename KeyType, typename MappableType1, typename MappableType2>
bool
has_contained(const MappableMap<KeyType, MappableType1>& mappables, const MappableType2& mappable)
{
    return std::any_of(std::cbegin(mappables), std::cend(mappables),
                       [&mappable] (const auto& p) {
                           return p.second.has_contained(mappable);
                       });
}

template <typename KeyType, typename MappableType1, typename MappableType2>
typename MappableFlatMultiSet<MappableType1>::size_type
count_contained(const MappableMap<KeyType, MappableType1>& mappables, const MappableType2& mappable)
{
    using SizeType = typename MappableFlatMultiSet<MappableType1>::size_type;
    
    return std::accumulate(std::cbegin(mappables), std::cend(mappables), SizeType {0},
                           [&mappable] (const auto curr, const auto& p) {
                               return curr + p.second.count_contained(mappable);
                           });
}

template <typename KeyType, typename MappableType1, typename MappableType2, typename MappableType3>
bool has_shared(const MappableMap<KeyType, MappableType1>& mappables,
                const MappableType2& mappable1, const MappableType3& mappable2)
{
    return std::any_of(std::cbegin(mappables), std::cend(mappables),
                       [&mappable1, &mappable2] (const auto& p) {
                           return p.second.has_shared(mappable1, mappable2);
                       });
}

template <typename KeyType, typename MappableType1, typename MappableType2, typename MappableType3>
typename MappableFlatMultiSet<MappableType1>::size_type
count_shared(const MappableMap<KeyType, MappableType1>& mappables, const MappableType2& mappable1,
             const MappableType3& mappable2)
{
    using SizeType = typename MappableFlatMultiSet<MappableType1>::size_type;
    
    return std::accumulate(std::cbegin(mappables), std::cend(mappables), SizeType {0},
                           [&mappable1, &mappable2] (const auto curr, const auto& p) {
                               return curr + p.second.count_shared(mappable1, mappable2);
                           });
}

template <typename ForwardIterator, typename KeyType, typename MappableType1, typename MappableType2>
ForwardIterator
find_first_shared(const MappableMap<KeyType, MappableType1>& mappables, ForwardIterator first,
                  ForwardIterator last, const MappableType2& mappable)
{
    using std::cbegin; using std::cend; using std::begin; using std::end;
    
    if (mappables.empty()) return last;
    
    if (mappables.size() == 1) {
        return find_first_shared(cbegin(mappables)->second, first, last, mappable);
    }
    
    std::vector<ForwardIterator> smallest(mappables.size());
    
    std::transform(cbegin(mappables), cend(mappables), begin(smallest),
                   [first, last, &mappable] (const auto& p) {
                       return find_first_shared(p.second, first, last, mappable);
                   });
    
    smallest.erase(std::remove_if(begin(smallest), end(smallest),
                                  [last] (auto it) { return it == last; }),
                   end(smallest));
    
    if (smallest.empty()) return last;
    
    return *std::min_element(cbegin(smallest), cend(smallest), []
                             (auto lhs, auto rhs) { return *lhs < *rhs; });
}

template <typename KeyType, typename MappableType, typename ForwardIterator>
size_t
max_count_if_shared_with_first(const MappableMap<KeyType, MappableType>& mappables,
                               ForwardIterator first, ForwardIterator last)
{
    size_t maximum {0};
    size_t count {0};
    
    for (const auto& map_pair : mappables) {
        count = count_if_shared_with_first(map_pair.second, first, last);
        if (count > maximum) {
            maximum = count;
        }
    }
    
    return maximum;
}

template <typename KeyType, typename MappableType1, typename MappableType2>
typename MappableFlatMultiSet<MappableType1>::const_iterator
leftmost_overlapped(const MappableMap<KeyType, MappableType1>& mappables, const MappableType2& mappable)
{
    if (mappables.empty()) {
        throw std::logic_error {"leftmost_overlapped called with empty MappableMap"};
    }
    
    auto first = std::cbegin(mappables);
    auto last  = std::cend(mappables);
    
    auto result = std::cbegin(first->second);
    
    while (first != last) {
        if (!first->second.empty()) {
            const auto overlapped = first->second.overlap_range(mappable);
            if (!overlapped.empty()) {
                result = std::cbegin(overlapped).base();
                ++first;
                break;
            }
        }
        ++first;
    }
    
    std::for_each(first, last, [&mappable, &result] (const auto& p) {
        const auto overlapped = p.second.overlap_range(mappable);
        if (!overlapped.empty() && begins_before(overlapped.front(), *result)) {
            result = std::cbegin(overlapped).base();
        }
    });
    
    return result;
}

template <typename KeyType, typename MappableType1, typename MappableType2>
typename MappableFlatMultiSet<MappableType1>::const_iterator
rightmost_overlapped(const MappableMap<KeyType, MappableType1>& mappables, const MappableType2& mappable)
{
    if (mappables.empty()) {
        throw std::logic_error {"rightmost_overlapped called with empty MappableMap"};
    }
    
    auto first = std::cbegin(mappables);
    auto last  = std::cend(mappables);
    
    auto result = std::cend(first->second);
    
    while (first != last) {
        if (!first->second.empty()) {
            const auto overlapped = first->second.overlap_range(mappable);
            if (!overlapped.empty()) {
                result = rightmost_mappable(overlapped).base();
                ++first;
                break;
            }
        }
        ++first;
    }
    
    std::for_each(first, last, [&mappable, &result] (const auto& p) {
        const auto overlapped = p.second.overlap_range(mappable);
        if (!overlapped.empty()) {
            result = rightmost_mappable(overlapped).base();
        }
    });
    
    return result;
}

template <typename KeyType, typename MappableType>
std::vector<unsigned>
positional_coverage(const MappableMap<KeyType, MappableType>& mappables, const GenomicRegion& region)
{
    std::vector<unsigned> result(region_size(region), 0);
    
    for (const auto& p : mappables) {
        const auto pcoverage = positional_coverage(p.second, region);
        
        std::transform(std::cbegin(result), std::cend(result), std::cbegin(pcoverage),
                       std::begin(result), std::plus<void> {});
    }
    
    return result;
}

template <typename KeyType, typename MappableType>
std::vector<unsigned>
positional_coverage(const MappableMap<KeyType, MappableType>& mappables)
{
    return positional_coverage(mappables, encompassing_region(mappables));
}

template <typename KeyType, typename MappableType1, typename MappableType2>
MappableMap<KeyType, MappableType1>
copy_overlapped(const MappableMap<KeyType, MappableType1>& mappables, const MappableType2& mappable)
{
    MappableMap<KeyType, MappableType1> result {mappables.size()};
    
    for (const auto& p : mappables) {
        result.emplace(p.first, copy_overlapped(p.second, mappable));
    }
    
    return result;
}

template <typename KeyType, typename MappableType1, typename MappableType2>
MappableMap<KeyType, MappableType1>
copy_contained(const MappableMap<KeyType, MappableType1>& mappables, const MappableType2& mappable)
{
    MappableMap<KeyType, MappableType1> result {mappables.size()};
    
    for (const auto& p : mappables) {
        result.emplace(p.first, copy_contained(p.second, mappable));
    }
    
    return result;
}

template <typename KeyType, typename MappableType1, typename MappableType2>
void erase_overlapped(MappableMap<KeyType, MappableType1>& mappables, const MappableType2& mappable)
{
    for (auto& p : mappables) {
        p.second.erase_overlapped(mappable);
    }
}

template <typename KeyType, typename MappableType1, typename MappableType2>
void erase_contained(MappableMap<KeyType, MappableType1>& mappables, const MappableType2& mappable)
{
    for (auto& p : mappables) {
        p.second.erase_contained(mappable);
    }
}

#endif
