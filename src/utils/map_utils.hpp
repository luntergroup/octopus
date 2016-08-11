// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef map_utils_hpp
#define map_utils_hpp

#include <algorithm>
#include <iterator>
#include <vector>
#include <utility>

#include "type_tricks.hpp"

namespace octopus {

template <typename MapType, typename T>
auto insert_or_assign(MapType& map, const typename MapType::key_type& k, T&& obj)
{
    const auto p = map.insert(std::make_pair(k, std::forward<T>(obj)));
    if (!p.second) {
        p.first->second = std::forward<T>(obj);
    }
    return p;
}

template <typename MapType>
auto extract_keys(const MapType& map)
{
    static_assert(is_map<MapType>, "MapType must be a map type");
    
    std::vector<typename MapType::key_type> result {};
    result.reserve(map.size());
    
    std::transform(std::cbegin(map), std::cend(map), std::back_inserter(result),
                   [] (const auto& p) { return p.first; });
    
    return result;
}

template <typename MapType>
auto extract_values(const MapType& map)
{
    static_assert(is_map<MapType>, "MapType must be a map type");
    
    std::vector<typename MapType::mapped_type> result {};
    result.reserve(map.size());
    
    std::transform(std::cbegin(map), std::cend(map), std::back_inserter(result),
                   [] (const auto& p) { return p.second; });
    
    return result;
}

template <typename MapType>
auto extract_value_sorted_keys(const MapType& map)
{
    static_assert(is_map<MapType>, "MapType must be a map type");
    
    std::vector<std::pair<typename MapType::mapped_type, typename MapType::key_type>> pairs {};
    pairs.reserve(map.size());
    
    std::transform(std::cbegin(map), std::cend(map), std::back_inserter(pairs),
                   [] (const auto& p) { return std::make_pair(p.second, p.first); });
    
    std::sort(std::begin(pairs), std::end(pairs),
              [] (const auto& lhs, const auto& rhs) { return lhs.first > rhs.first; });
    
    std::vector<typename MapType::key_type> result {};
    result.reserve(pairs.size());
    
    std::transform(std::cbegin(pairs), std::cend(pairs), std::back_inserter(result),
                   [] (const auto& p) { return p.second; });
    
    return result;
}

template <typename MapType, typename UnaryFunction>
void for_each_value(MapType& map, UnaryFunction f)
{
    static_assert(is_map<MapType>, "map must be a map type");
    
    for (auto& p : map) {
        f(*p.second);
    }
}

} // namespace octopus

#endif
