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

namespace detail {

template <typename T, typename = void>
struct IsIndexedHelper : std::false_type {};

template <typename T>
struct IsIndexedHelper<T,
        std::enable_if_t<std::is_integral<decltype(index_of(std::declval<T>()))>::value>
    > : std::true_type {};

} // namespace detail

template <typename T>
constexpr bool is_indexed = detail::IsIndexedHelper<T>::value;

} // namespace octopus

#endif
