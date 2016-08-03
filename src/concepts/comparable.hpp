//
//  comparable.hpp
//  Octopus
//
//  Created by Daniel Cooke on 13/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_comparable_hpp
#define Octopus_comparable_hpp

namespace octopus {

/**
 A class that is derived from Comparable must implement operator== and operator< and will 
 then have all other comparison operators defined.
 */
template <typename T>
class Comparable {};

template <typename T>
inline bool operator!=(const Comparable<T>& lhs, const Comparable<T>& rhs)
{
    return !operator==(static_cast<const T&>(lhs), static_cast<const T&>(rhs));
}

template <typename T>
inline bool operator>(const Comparable<T>& lhs, const Comparable<T>& rhs)
{
    return operator<(static_cast<const T&>(rhs), static_cast<const T&>(lhs));
}

template <typename T>
inline bool operator<=(const Comparable<T>& lhs, const Comparable<T>& rhs)
{
    return !operator>(static_cast<const T&>(lhs), static_cast<const T&>(rhs));
}

template <typename T>
inline bool operator>=(const Comparable<T>& lhs, const Comparable<T>& rhs)
{
    return !operator<(static_cast<const T&>(lhs), static_cast<const T&>(rhs));
}

} // namespace octopus

#endif
