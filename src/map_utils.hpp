//
//  map_utils.hpp
//  Octopus
//
//  Created by Daniel Cooke on 28/04/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_map_utils_hpp
#define Octopus_map_utils_hpp

#include <cstddef>   // std::size_t
#include <algorithm> // std::any_of, std::transform, std::min_element
#include <iterator>  // std::distance, std::cbegin, std::cend, std::prev
#include <vector>

#include "genomic_region.hpp"
#include "mappable.hpp"
#include "mappable_algorithms.hpp"

template <typename MappableMap, typename MappableType>
std::size_t count_overlapped(const MappableMap& mappables, const MappableType& m,
                             MappableRangeOrder order=MappableRangeOrder::ForwardSorted)
{
    std::size_t result {};
    
    for (const auto& map_pair : mappables) {
        result += count_overlapped(std::cbegin(map_pair.second), std::cend(map_pair.second), m, order);
    }
    
    return result;
}

template <typename MappableMap, typename MappableType1, typename MappableType2>
bool has_shared(const MappableMap& mappables, const MappableType1& lhs, const MappableType2& rhs,
                MappableRangeOrder order=MappableRangeOrder::ForwardSorted)
{
    return std::any_of(std::cbegin(mappables), std::cend(mappables),
                       [&lhs, &rhs, order] (const auto& map_pair) {
                           return has_shared(std::cbegin(map_pair.second), std::cend(map_pair.second),
                                             lhs, rhs, order);
                       });
}

template <typename MappableMap, typename MappableType1, typename MappableType2>
std::size_t
count_shared(const MappableMap& mappables, const MappableType1& lhs, const MappableType2& rhs,
             MappableRangeOrder order=MappableRangeOrder::ForwardSorted)
{
    std::size_t result {};
    
    for (const auto& map_pair : mappables) {
        result += count_shared(std::cbegin(map_pair.second), std::cend(map_pair.second), lhs, rhs, order);
    }
    
    return result;
}

/**
 Works the same way as the other find_first_shared, but operators over a map of ranges
 (e.g. std::map or std::unordered_map).
 */
template <typename ForwardIterator, typename MappableMap, typename MappableType>
ForwardIterator
find_first_shared(const MappableMap& mappables, ForwardIterator first, ForwardIterator last,
                  const MappableType& mappable, MappableRangeOrder order = MappableRangeOrder::ForwardSorted)
{
    if (mappables.empty()) return last;
    
    if (mappables.size() == 1) {
        return find_first_shared(std::cbegin(std::cbegin(mappables)->second),
                                 std::cend(std::cbegin(mappables)->second), first, last, mappable, order);
    }
    
    std::vector<ForwardIterator> smallest {};
    smallest.reserve(mappables.size());
    
    std::transform(std::cbegin(mappables), std::cend(mappables), std::back_inserter(smallest),
                   [first, last, &mappable, order] (const auto& p) {
                       return find_first_shared(std::cbegin(p.second), std::cend(p.second), first,
                                                last, mappable, order);
                   });
    
    return *std::min_element(std::cbegin(smallest), std::cend(smallest), []
                             (ForwardIterator lhs, ForwardIterator rhs) {
                                 return *lhs < *rhs;
                             });
}

template <typename MappableMap, typename ForwardIterator>
std::size_t
max_count_if_shared_with_first(const MappableMap& mappables, ForwardIterator first, ForwardIterator last,
                               MappableRangeOrder order=MappableRangeOrder::ForwardSorted)
{
    std::size_t maximum {0};
    std::size_t count {};
    
    for (const auto& map_pair : mappables) {
        count = count_if_shared_with_first(std::cbegin(map_pair.second), std::cend(map_pair.second),
                                           first, last, order);
        if (count > maximum) {
            maximum = count;
        }
    }
    
    return maximum;
}

template <typename MappableMap>
typename MappableMap::mapped_type::const_iterator
leftmost_sorted_mappable(const MappableMap& mappables)
{
    using Iterator = typename MappableMap::mapped_type::const_iterator;
    
    Iterator result {std::cbegin(std::cbegin(mappables)->second)};
    
    for (const auto& map_pair : mappables) {
        if (!map_pair.second.empty() && begins_before(*std::cbegin(map_pair.second), *result)) {
            result = std::cbegin(map_pair.second);
        }
    }
    
    return result;
}

template <typename MappableMap, typename MappableType>
typename MappableMap::mapped_type::const_iterator
leftmost_overlapped(const MappableMap& mappables, const MappableType& mappable,
                    MappableRangeOrder order=MappableRangeOrder::ForwardSorted)
{
    using Iterator = typename MappableMap::mapped_type::const_iterator;
    
    Iterator result;
    
    // To find a default value
    for (const auto& map_pair : mappables) {
        if (map_pair.second.size() > 0) {
            result = std::prev(std::cend(map_pair.second));
            break;
        }
    }
    
    for (const auto& map_pair : mappables) {
        auto overlapped = overlap_range(std::cbegin(map_pair.second), std::cend(map_pair.second), mappable, order);
        
        if (!overlapped.empty() && overlapped.begin() != std::cend(map_pair.second) && begins_before(overlapped.front(), *result)) {
            result = overlapped.begin().base();
        }
    }
    
    return result;
}

template <typename MappableMap, typename MappableType>
typename MappableMap::mapped_type::const_iterator
rightmost_overlapped(const MappableMap& mappables, const MappableType& mappable,
                     MappableRangeOrder order=MappableRangeOrder::ForwardSorted)
{
    using Iterator = typename MappableMap::mapped_type::const_iterator;
    
    Iterator result {std::cbegin(std::cbegin(mappables)->second)};
    
    for (const auto& map_pair : mappables) {
        auto overlapped = overlap_range(std::cbegin(map_pair.second), std::cend(map_pair.second), mappable, order);
        
        if (!overlapped.empty() && overlapped.begin() != overlapped.end() && ends_before(*result, overlapped.back())) {
            result = std::prev(overlapped.end().base());
        }
    }
    
    return result;
}

template <typename MapType>
std::vector<typename MapType::key_type>
value_sorted_keys(const MapType& map)
{
    std::vector<std::pair<typename MapType::mapped_type, typename MapType::key_type>> pairs {};
    pairs.reserve(map.size());
    
    std::transform(std::cbegin(map), std::cend(map), std::back_inserter(pairs),
                   [] (const auto& p) { return std::make_pair(p.second, p.first); });
    
    std::sort(pairs.begin(), pairs.end(), [] (const auto& lhs, const auto& rhs) { return lhs.first > rhs.first; });
    
    std::vector<typename MapType::key_type> result {};
    result.reserve(pairs.size());
    
    std::transform(pairs.cbegin(), pairs.cend(), std::back_inserter(result),
                   [] (const auto& p) { return p.second; });
    
    return result;
}

#endif
