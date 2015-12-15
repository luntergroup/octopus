//
//  mappable_debug.hpp
//  Octopus
//
//  Created by Daniel Cooke on 14/12/2015.
//  Copyright © 2015 Oxford University. All rights reserved.
//

#include <iostream>
#include <algorithm>
#include <iterator>

#include "mappable_algorithms.hpp"
#include "mappable_set.hpp"
#include "mappable_map.hpp"

#ifndef mappable_debug_hpp
#define mappable_debug_hpp

template <typename T, typename MappableType>
void print_overlapped(const MappableSet<T>& mappables, const MappableType& mappable,
                      const std::string& delim = "\n")
{
    const auto overlapped = mappables.overlap_range(mappable);
    std::copy(std::cbegin(overlapped), std::cend(overlapped),
              std::ostream_iterator<T>(std::cout, delim.c_str()));
}

template <typename Container, typename MappableType>
void print_overlapped(const Container& mappables, const MappableType& mappable,
                      const std::string& delim = "\n")
{
    using T = typename std::iterator_traits<typename Container::const_iterator>::value_type;
    const auto overlapped = overlap_range(mappables, mappable);
    std::copy(std::cbegin(overlapped), std::cend(overlapped),
              std::ostream_iterator<T>(std::cout, delim.c_str()));
}

template <typename K, typename T, typename MappableType>
void print_overlapped(const MappableMap<K, T>& mappables, const MappableType& mappable,
                      const std::string& delim = "\n")
{
    for (const auto& p : mappables) {
        std::cout << p.first << std::endl;
        print_overlapped(p.second, mappable, delim);
    }
}

template <typename T, typename MappableType>
void print_contained(const MappableSet<T>& mappables, const MappableType& mappable,
                     const std::string& delim = "\n")
{
    const auto contained = mappables.contained_range(mappable);
    std::copy(std::cbegin(contained), std::cend(contained),
              std::ostream_iterator<T>(std::cout, delim.c_str()));
}

template <typename Container, typename MappableType>
void print_contained(const Container& mappables, const MappableType& mappable,
                     const std::string& delim = "\n")
{
    using T = typename std::iterator_traits<typename Container::const_iterator>::value_type;
    const auto contained = contained_range(mappables, mappable);
    std::copy(std::cbegin(contained), std::cend(contained),
              std::ostream_iterator<T>(std::cout, delim.c_str()));
}

template <typename K, typename T, typename MappableType>
void print_contained(const MappableMap<K, T>& mappables, const MappableType& mappable,
                     const std::string& delim = "\n")
{
    for (const auto& p : mappables) {
        std::cout << p.first << std::endl;
        print_contained(p.second, mappable, delim);
    }
}

#endif /* mappable_debug_hpp */
