//
//  map_utils.h
//  Octopus
//
//  Created by Daniel Cooke on 28/04/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_map_utils_h
#define Octopus_map_utils_h

#include <cstddef>   // std::size_t
#include <algorithm> // std::any_of, std::transform, std::min_element
#include <iterator>  // std::distance, std::cbegin, std::cend, std::prev
#include <vector>

#include "genomic_region.h"
#include "mappable.h"

template <typename MapType>
inline
typename MapType::mapped_type
sum(const MapType& map)
{
    typename MapType::mapped_type result {};
    
    for (const auto& map_pair : map) {
        result += map_pair.second;
    }
    
    return result;
}

template <typename ResultType, typename MapType, typename UnaryOperation>
inline
ResultType
sum(const MapType& map, UnaryOperation op)
{
    ResultType result {};
    
    for (const auto& map_pair : map) {
        result += op(map_pair.second);
    }
    
    return result;
}

template <typename MappableMap, typename MappableType>
inline
std::size_t count_overlapped(const MappableMap& mappables, const MappableType& m)
{
    std::size_t result {};
    
    for (const auto& map_pair : mappables) {
        result += count_overlapped(std::cbegin(map_pair.second), std::cend(map_pair.second), m);
    }
    
    return result;
}

template <typename MappableMap, typename MappableType1, typename MappableType2>
inline
bool has_shared(const MappableMap& mappables, const MappableType1& lhs, const MappableType2& rhs)
{
    return std::any_of(std::cbegin(mappables), std::cend(mappables),
                       [&lhs, &rhs] (const auto& map_pair) {
                           return has_shared(std::cbegin(map_pair.second), std::cend(map_pair.second),
                                             lhs, rhs);
                       });
}

template <typename MappableMap, typename MappableType1, typename MappableType2>
inline
std::size_t
count_shared(const MappableMap& mappables, const MappableType1& lhs, const MappableType2& rhs)
{
    std::size_t result {};
    
    for (const auto& map_pair : mappables) {
        result += count_shared(std::cbegin(map_pair.second), std::cend(map_pair.second), lhs, rhs);
    }
    
    return result;
}

/**
 Works the same way as the other find_first_shared, but operators over a map of ranges
 (e.g. std::map or std::unordered_map).
 */
template <typename ForwardIterator, typename MappableMap, typename MappableType>
ForwardIterator find_first_shared(const MappableMap& mappables, ForwardIterator first,
                                  ForwardIterator last, const MappableType& mappable)
{
    if (mappables.empty()) return last;
    
    if (mappables.size() == 1) {
        return find_first_shared(std::cbegin(std::cbegin(mappables)->second),
                                 std::cend(std::cbegin(mappables)->second), first, last, mappable);
    }
    
    std::vector<ForwardIterator> smallest {};
    smallest.reserve(mappables.size());
    
    std::transform(std::cbegin(mappables), std::cend(mappables), std::back_inserter(smallest),
                   [first, last, &mappable] (const auto& p) {
                       return find_first_shared(std::cbegin(p.second), std::cend(p.second), first,
                                                last, mappable);
                   });
    
    return *std::min_element(std::cbegin(smallest), std::cend(smallest), []
                             (ForwardIterator lhs, ForwardIterator rhs) {
                                 return *lhs < *rhs;
                             });
}

template <typename MappableMap, typename ForwardIterator>
inline
std::size_t max_count_if_shared_with_first(const MappableMap& mappables,
                                           ForwardIterator first, ForwardIterator last)
{
    std::size_t maximum {0};
    std::size_t count {};
    
    for (const auto& map_pair : mappables) {
        count = count_if_shared_with_first(std::cbegin(map_pair.second), std::cend(map_pair.second),
                                           first, last);
        if (count > maximum) {
            maximum = count;
        }
    }
    
    return maximum;
}

template <typename MappableMap>
inline
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
inline
typename MappableMap::mapped_type::const_iterator
leftmost_overlapped(const MappableMap& mappables, const MappableType& mappable)
{
    using Iterator = typename MappableMap::mapped_type::const_iterator;
    
    Iterator result;
    std::pair<Iterator, Iterator> overlapped_mappables;
    
    // To find a default value
    for (const auto& map_pair : mappables) {
        if (map_pair.second.size() > 0) {
            result = std::prev(std::cend(map_pair.second));
            break;
        }
    }
    
    for (const auto& map_pair : mappables) {
        overlapped_mappables = overlap_range(std::cbegin(map_pair.second),
                                             std::cend(map_pair.second), mappable);
        
        if (overlapped_mappables.first != std::cend(map_pair.second) &&
            begins_before(*overlapped_mappables.first, *result)) {
            result = overlapped_mappables.first;
        }
    }
    
    return result;
}

template <typename MappableMap, typename MappableType>
inline
typename MappableMap::mapped_type::const_iterator
rightmost_overlapped(const MappableMap& mappables, const MappableType& mappable)
{
    using Iterator = typename MappableMap::mapped_type::const_iterator;
    
    Iterator result {std::cbegin(std::cbegin(mappables)->second)};
    std::pair<Iterator, Iterator> overlapped_mappables;
    
    for (const auto& map_pair : mappables) {
        overlapped_mappables = overlap_range(std::cbegin(map_pair.second),
                                             std::cend(map_pair.second), mappable);
        
        if (overlapped_mappables.first != overlapped_mappables.second &&
            ends_before(*result, *std::prev(overlapped_mappables.second))) {
            result = std::prev(overlapped_mappables.second);
        }
    }
    
    return result;
}

#endif
