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

#endif /* defined(__Octopus__utils__) */
