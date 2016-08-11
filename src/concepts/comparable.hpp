// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef comparable_hpp
#define comparable_hpp

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
