// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef array_tricks_hpp
#define array_tricks_hpp

#include <array>
#include <type_traits>
#include <utility>

namespace octopus {

namespace detail {
template <typename T, std::size_t...Is>
constexpr std::array<T, sizeof...(Is)>
make_array(const T& value, std::index_sequence<Is...>)
{
    return {{(static_cast<void>(Is), value)...}};
}
} // namespace detail

template <std::size_t N, typename T>
constexpr std::array<T, N> make_array(const T& value)
{
    return detail::make_array(value, std::make_index_sequence<N>());
}

namespace detail {
template <typename T, std::size_t...Is>
constexpr std::array<T, sizeof...(Is) + 1>
make_array(const T& first, const T& rest, std::index_sequence<Is...>)
{
    return {{first, (static_cast<void>(Is), rest)...}};
}
} // namespace detail

template <std::size_t N, typename T>
constexpr std::array<T, N> make_array(const T& first, const T& rest)
{
    return detail::make_array(first, rest, std::make_index_sequence<N - 1>());
}
template <typename T>
constexpr std::array<T, 1> make_array(const T& first, const T& rest)
{
    return make_array<1>(first);
}

namespace detail {
template <typename UnaryOperation, typename T, std::size_t N, std::size_t...Is>
constexpr auto
transform(UnaryOperation&& op, const std::array<T, N>& values, std::index_sequence<Is...>)
-> std::array<decltype(op(T{})), N>
{
    return {{op(std::get<Is>(values))...}};
}
} // namespace detail
template <typename T, std::size_t N, typename UnaryOperation>
constexpr auto
transform(UnaryOperation&& op, const std::array<T, N>& values)
{
    return detail::transform(std::forward<UnaryOperation>(op), values, std::make_index_sequence<N>());
}

namespace detail {
template <typename BinaryOperation, typename T, std::size_t N, std::size_t...Is>
constexpr auto
transform(BinaryOperation&& op, const std::array<T, N>& lhs, const std::array<T, N>& rhs, std::index_sequence<Is...>)
-> std::array<decltype(op(T{}, T{})), N>
{
    return {{op(std::get<Is>(lhs), std::get<Is>(rhs))...}};
}
} // namespace detail
template <typename T, std::size_t N, typename BinaryOperation>
constexpr auto
transform(BinaryOperation&& op, const std::array<T, N>& lhs, const std::array<T, N>& rhs)
{
    return detail::transform(std::forward<BinaryOperation>(op), lhs, rhs, std::make_index_sequence<N>());
}

namespace detail {
template <typename BinaryOperation, typename T, std::size_t N, std::size_t...Is>
constexpr void
adjacent_apply(BinaryOperation&& op, const std::array<T, N>& source, std::array<T, N>& dest, std::index_sequence<Is...>)
{
    (void) std::initializer_list<int>{(std::get<Is>(dest) = op(std::get<Is>(source), std::get<Is + 1>(source)), 0)...};
}
} // namespace detail
template <typename BinaryOperation, typename T, std::size_t N>
constexpr void
adjacent_apply(BinaryOperation&& op, const std::array<T, N>& source,  std::array<T, N>& dest)
{
    return detail::adjacent_apply(std::forward<BinaryOperation>(op), source, dest, std::make_index_sequence<N - 1>());
}

namespace detail {
template <typename BinaryOperation, typename T, std::size_t N, std::size_t...Is>
constexpr void
adjacent_apply_reverse(BinaryOperation&& op, const std::array<T, N>& source, std::array<T, N>& dest,
                       std::index_sequence<Is...>)
{
    (void) std::initializer_list<int>{
    (std::get<N - Is - 1>(dest) = op(std::get<N - Is - 2>(source), std::get<N - Is - 1>(source)), 0)...};
}
} // namespace detail
template <typename BinaryOperation, typename T, std::size_t N>
constexpr void
adjacent_apply_reverse(BinaryOperation&& op, const std::array<T, N>& source,  std::array<T, N>& dest)
{
    return detail::adjacent_apply_reverse(std::forward<BinaryOperation>(op), source, dest, std::make_index_sequence<N - 1>());
}

} // namespace octopus

#endif
