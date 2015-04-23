//
//  test_utils.h
//  Octopus
//
//  Created by Daniel Cooke on 23/04/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_test_utils_h
#define Octopus_test_utils_h

#include <cmath>
#include <string>
#include <unordered_map>

#include "haplotype.h"

template <typename RealType>
bool is_close_to_one(RealType x)
{
    return std::abs(x - 1) < 0.0000000000001;
}

inline std::string get_label(const std::unordered_map<Haplotype, std::string>& labels, const Haplotype& haplotype)
{
    return (labels.count(haplotype) > 0) ? labels.at(haplotype) : "other";
}

#endif
