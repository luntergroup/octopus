//
//  equitable.hpp
//  Octopus
//
//  Created by Daniel Cooke on 13/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_equitable_h
#define Octopus_equitable_h

namespace octopus {

template <typename T>
class Equitable {};

template <typename T>
inline bool operator!=(const Equitable<T>& lhs, const Equitable<T>& rhs)
{
    return !operator==(static_cast<const T&>(lhs), static_cast<const T&>(rhs));
}

} // namespace octopus

#endif
