//
//  maths.h
//  Octopus
//
//  Created by Daniel Cooke on 23/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__maths__
#define __Octopus__maths__

#include <vector>
#include <numeric>

template <typename T>
T multinomial_coefficient(const std::vector<T>& xs)
{
    auto x_0 = std::accumulate(xs.cbegin(), xs.cend(), T {});
    return x_0;
}

#endif /* defined(__Octopus__maths__) */
