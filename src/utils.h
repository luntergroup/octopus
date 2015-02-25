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
#include <sstream>
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
