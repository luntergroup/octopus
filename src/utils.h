//
//  utils.h
//  Octopus
//
//  Created by Daniel Cooke on 10/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__utils__
#define __Octopus__utils__

#include <vector>
#include <string>
#include <cstring>   // std::strlen
#include <sstream>
#include <algorithm> // std::equal, std::swap
#include <iterator>  // std::next
#include <cstddef>   // std::size_t

template <typename T>
std::vector<std::string> split(T&& s, char delim) {
    std::vector<std::string> elems;
    std::stringstream ss(s);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.emplace_back(item);
    }
    return elems;
}

template <typename T>
bool is_prefix(T this_sequence, T that_sequence)
{
    return std::equal(std::cbegin(this_sequence), std::cend(this_sequence),
                      std::cbegin(that_sequence));
}

template <typename T>
bool is_suffix(T this_sequence, T that_sequence)
{
    return std::equal(std::cbegin(this_sequence), std::cend(this_sequence),
                      std::next(std::cbegin(that_sequence)));
}

inline std::size_t stringlen(const char* str)
{
    return std::strlen(str);
}

inline std::size_t stringlen(const std::string& str)
{
    return str.size();
}

template <typename C>
inline C kth_permutation(std::size_t k, C container)
{
    for (std::size_t j = 1; j < container.size(); ++j) {
        std::swap(container[k % (j + 1)], container[j]);
        k = k / (j + 1);
    }
    return container;
}

#endif /* defined(__Octopus__utils__) */
