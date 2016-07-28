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
#include <cstring>
#include <sstream>
#include <algorithm>
#include <iterator>
#include <cstddef>
#include <type_traits>
#include <iomanip>

namespace octopus
{
template <typename T>
std::vector<std::string> split(const T& str, const char delim) {
    std::vector<std::string> elems;
    elems.reserve(std::count(std::cbegin(str), std::cend(str), delim) + 1);
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

inline bool find(const std::string& lhs, const std::string& rhs)
{
    return lhs.find(rhs) != std::string::npos;
}

template <typename T, typename = typename std::enable_if_t<std::is_floating_point<T>::value>>
std::string to_string(const T val, const unsigned precision = 2)
{
    std::ostringstream out;
    out << std::fixed << std::setprecision(precision) << val;
    return out.str();
}

template <typename T>
std::vector<std::string> to_strings(const std::vector<T>& values)
{
    std::vector<std::string> result {};
    result.reserve(values.size());
    std::transform(std::cbegin(values), std::cend(values), std::back_inserter(result),
                   [] (auto value) { return std::to_string(value); });
    return result;
}

template <typename T, typename = typename std::enable_if_t<std::is_floating_point<T>::value>>
std::vector<std::string> to_strings(const std::vector<T>& values, const unsigned precision = 2)
{
    std::vector<std::string> result {};
    result.reserve(values.size());
    std::transform(std::cbegin(values), std::cend(values), std::back_inserter(result),
                   [precision] (auto value) { return to_string(value, precision); });
    return result;
}

} // namespace octopus

#endif /* defined(__Octopus__string_utils__) */
