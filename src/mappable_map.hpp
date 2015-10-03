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
#include <algorithm> // std::for_each, std::any_of, std::transform
#include <numeric>   // std::accumulate
#include <iterator>  // std::begin, std::end, std::make_move_iterator
#include <stdexcept>
#include <cstddef>   // size_t

#include "mappable_set.hpp"

template <typename KeyType, typename MappableType>
using MappableMap = std::unordered_map<KeyType, MappableSet<MappableType>>;

template <typename Map>
MappableMap<typename Map::key_type, typename Map::mapped_type::value_type>
make_mappable_map(Map map)
{
    using T = typename Map::mapped_type::value_type;
    
    MappableMap<typename Map::key_type, T> result {};
    result.reserve(map.size());
    
    using std::make_move_iterator; using std::begin; using std::end;
    
    std::for_each(make_move_iterator(begin(map)), make_move_iterator(end(map)), [&result] (auto&& p) {
        result.emplace(std::move(p.first),
                       MappableSet<T> { make_move_iterator(begin(p.second)), make_move_iterator(end(p.second)) });
    });
    
    return result;
}

template <typename KeyType, typename MappableType>
size_t count_mappables(const MappableMap<KeyType, MappableType>& map)
{
    return std::accumulate(std::cbegin(map), std::cend(map), size_t {},
                           [] (size_t prev, const auto& v) { return prev + v.second.size(); });
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
typename MappableSet<MappableType1>::size_type
count_overlapped(const MappableMap<KeyType, MappableType1>& mappables, const MappableType2& mappable)
{
    using SizeType = typename MappableSet<MappableType1>::size_type;
    return std::accumulate(std::cbegin(mappables), std::cend(mappables), SizeType {},
                           [&mappable] (SizeType x, const auto& p) {
                               return x + p.second.count_overlapped(mappable);
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
typename MappableSet<MappableType1>::size_type
count_contained(const MappableMap<KeyType, MappableType1>& mappables, const MappableType2& mappable)
{
    using SizeType = typename MappableSet<MappableType1>::size_type;
    return std::accumulate(std::cbegin(mappables), std::cend(mappables), SizeType {},
                           [&mappable] (SizeType x, const auto& p) {
                               return x + p.second.count_contained(mappable);
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
typename MappableSet<MappableType1>::size_type
count_shared(const MappableMap<KeyType, MappableType1>& mappables, const MappableType2& mappable1,
             const MappableType3& mappable2)
{
    using SizeType = typename MappableSet<MappableType1>::size_type;
    return std::accumulate(std::cbegin(mappables), std::cend(mappables), SizeType {},
                           [&mappable1, &mappable2] (SizeType x, const auto& p) {
                               return x + p.second.count_shared(mappable1, mappable2);
                           });
}

template <typename ForwardIterator, typename KeyType, typename MappableType1, typename MappableType2>
ForwardIterator
find_first_shared(const MappableMap<KeyType, MappableType1>& mappables, ForwardIterator first,
                  ForwardIterator last, const MappableType2& mappable)
{
    if (mappables.empty()) return last;
    
    if (mappables.size() == 1) {
        return find_first_shared(mappables, first, last, mappable);
    }
    
    std::vector<ForwardIterator> smallest(mappables.size());
    
    std::transform(std::cbegin(mappables), std::cend(mappables), smallest.begin(),
                   [first, last, &mappables, &mappable] (const auto& p) {
                       return find_first_shared(mappables, first, last, mappable);
                   });
    
    return *std::min_element(std::cbegin(smallest), std::cend(smallest), []
                             (auto lhs, auto rhs) {
                                 return *lhs < *rhs;
                             });
}

template <typename KeyType, typename MappableType1, typename MappableType2>
typename MappableSet<MappableType1>::const_iterator
leftmost_overlapped(const MappableMap<KeyType, MappableType1>& mappables, const MappableType2& mappable)
{
    if (mappables.empty()) {
        throw std::runtime_error {"cannot find leftmost_overlapped of empty MappableMap"};
    }
    
    if (mappables.size() == 1) {
        return mappables.cbegin()->second.overlap_range(mappable).begin().base();
    }
    
    std::vector<typename MappableSet<MappableType1>::const_iterator> smallest(mappables.size());
    
    std::transform(std::cbegin(mappables), std::cend(mappables), smallest.begin(),
                   [&mappable] (const auto& p) {
                       return p.second.overlap_range(mappable).begin().base();
                   });
    
    return *std::min_element(std::cbegin(smallest), std::cend(smallest), []
                             (auto lhs, auto rhs) {
                                 return *lhs < *rhs;
                             });
}

template <typename KeyType, typename MappableType1, typename MappableType2>
typename MappableSet<MappableType1>::const_iterator
rightmost_overlapped(const MappableMap<KeyType, MappableType1>& mappables, const MappableType2& mappable)
{
    if (mappables.empty()) {
        throw std::runtime_error {"cannot find leftmost_overlapped of empty MappableMap"};
    }
    
    if (mappables.size() == 1) {
        auto overlapped = mappables.cbegin()->second.overlap_range(mappable);
        return rightmost_mappable(overlapped.begin(), overlapped.end()).base();
    }
    
    std::vector<typename MappableSet<MappableType1>::const_iterator> largest(mappables.size());
    
    std::transform(std::cbegin(mappables), std::cend(mappables), largest.begin(),
                   [&mappable] (const auto& p) {
                       auto overlapped = p.second.overlap_range(mappable);
                       return rightmost_mappable(overlapped.begin(), overlapped.end()).base();
                   });
    
    return *std::max_element(std::cbegin(largest), std::cend(largest), []
                             (auto lhs, auto rhs) {
                                 return ends_before(*lhs, *rhs);
                             });
}

#endif
