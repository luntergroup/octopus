// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef Octopus_equitable_h
#define Octopus_equitable_h

namespace octopus {

/**
 A class that is derived from Equitable must implement operator== and will then have
 operator!= defined.
 */
template <typename T>
class Equitable {};

template <typename T>
inline bool operator!=(const Equitable<T>& lhs, const Equitable<T>& rhs)
{
    return !operator==(static_cast<const T&>(lhs), static_cast<const T&>(rhs));
}

} // namespace octopus

#endif
