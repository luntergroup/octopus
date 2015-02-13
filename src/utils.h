//
//  utils.h
//  Octopus
//
//  Created by Daniel Cooke on 10/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__utils__
#define __Octopus__utils__

template <typename T>
class Comparable
{
};

template <typename T>
inline bool operator == (const Comparable<T>& lhs, const Comparable<T>& rhs)
{
    const T& d1 = static_cast<T const&>(lhs);
    const T& d2 = static_cast<T const&>(rhs);
    return !(d1 < d2) && !(d2 < d1);
}

#endif /* defined(__Octopus__utils__) */
