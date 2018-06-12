// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef equitable_hpp
#define equitable_hpp

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
