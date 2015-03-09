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
#include <algorithm> // std::equal
#include <iterator>  // std::next
#include <cstddef>   // std::size_t
#include <boost/utility/string_ref.hpp>
#include <boost/functional/hash.hpp>

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

namespace std
{
    template<>
    struct hash<boost::string_ref> {
        size_t operator()(boost::string_ref const& sr) const {
            return boost::hash_range(sr.begin(), sr.end());
        }
    };
}

#endif /* defined(__Octopus__utils__) */
