//
//  map_utils.hpp
//  Octopus
//
//  Created by Daniel Cooke on 28/04/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_map_utils_hpp
#define Octopus_map_utils_hpp

#include <algorithm>
#include <iterator>
#include <vector>

template <typename MapType>
std::vector<typename MapType::key_type> get_keys(const MapType& map)
{
    std::vector<typename MapType::key_type> result {};
    result.reserve(map.size());
    
    std::transform(std::cbegin(map), std::cend(map), std::back_inserter(result),
                   [] (const auto& p) { return p.first; });
    
    return result;
}

template <typename MapType>
std::vector<typename MapType::mapped_type> get_values(const MapType& map)
{
    std::vector<typename MapType::mapped_type> result {};
    result.reserve(map.size());
    
    std::transform(std::cbegin(map), std::cend(map), std::back_inserter(result),
                   [] (const auto& p) { return p.second; });
    
    return result;
}

template <typename MapType>
std::vector<typename MapType::key_type> get_value_sorted_keys(const MapType& map)
{
    std::vector<std::pair<typename MapType::mapped_type, typename MapType::key_type>> pairs {};
    pairs.reserve(map.size());
    
    std::transform(std::cbegin(map), std::cend(map), std::back_inserter(pairs),
                   [] (const auto& p) { return std::make_pair(p.second, p.first); });
    
    std::sort(pairs.begin(), pairs.end(),
              [] (const auto& lhs, const auto& rhs) { return lhs.first > rhs.first; });
    
    std::vector<typename MapType::key_type> result {};
    result.reserve(pairs.size());
    
    std::transform(std::cbegin(pairs), std::cend(pairs), std::back_inserter(result),
                   [] (const auto& p) { return p.second; });
    
    return result;
}

#endif
