//
//  string_utils.h
//  Octopus
//
//  Created by Daniel Cooke on 10/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__string_utils__
#define __Octopus__string_utils__

#include <vector>
#include <string>
#include <cstring>   // std::strlen
#include <sstream>   // std::stringstream
#include <algorithm> // std::equal, std::swap
#include <iterator>  // std::next
#include <cstddef>   // std::size_t

template <typename T>
std::vector<std::string> split(T&& str, char delim) {
    std::vector<std::string> elems;
    std::stringstream ss(str);
    std::string item;
    while (std::getline(ss, item, delim)) {
        elems.emplace_back(item);
    }
    return elems;
}

template <typename T>
bool is_prefix(const T& lhs, const T& rhs)
{
    return std::equal(std::cbegin(lhs), std::cend(lhs), std::cbegin(rhs));
}

template <typename T>
bool is_suffix(const T& lhs, const T& rhs)
{
    return std::equal(std::cbegin(lhs), std::cend(lhs), std::next(std::cbegin(rhs)));
}

inline std::size_t stringlen(const char* str)
{
    return std::strlen(str);
}

inline std::size_t stringlen(const std::string& str)
{
    return str.size();
}

inline bool contains(const std::string& lhs, const std::string& rhs)
{
    return lhs.find(rhs) != std::string::npos;
}

#endif /* defined(__Octopus__string_utils__) */
