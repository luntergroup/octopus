// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef indexed_hpp
#define indexed_hpp

#include <cstddef>
#include <type_traits>

namespace octopus {

template <typename T, typename IndexType = std::size_t>
class Indexed {};

template <typename T, typename IndexType>
IndexType index_of(const Indexed<T, IndexType>& indexed) noexcept
{
    return static_cast<const T&>(indexed).index();
}

template <typename T, typename = void>
struct is_indexed : std::false_type {};
template <typename T>
struct is_indexed<T, std::enable_if_t<std::is_integral<decltype(index_of(std::declval<T>()))>::value> > : std::true_type {};

template <typename T>
constexpr bool is_indexed_v = is_indexed<T>::value;

} // namespace octopus

#endif
