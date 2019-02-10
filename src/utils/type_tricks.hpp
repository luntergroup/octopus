// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef type_tricks_hpp
#define type_tricks_hpp

#include <type_traits>
#include <utility>

namespace octopus {

namespace detail {
    template <typename Container, typename = void>
    struct HasReserve : std::false_type {};
    
    template <typename Container>
    struct HasReserve<Container, std::enable_if_t<
            std::is_same<decltype(std::declval<Container>().reserve(std::declval<typename Container::size_type>())),
        void>::value>>
    : std::true_type
    {};
    
    template <typename Container>
    void reserve_if_enabled(Container& c, typename Container::size_type n, std::true_type)
    {
        c.reserve(n);
    }
    
    template <typename Container>
    void reserve_if_enabled(Container& c, typename Container::size_type n, std::false_type) {}
    
    template <typename Container, typename = void>
    struct HasShrinkToFit : std::false_type {};
    
    template <typename Container>
    struct HasShrinkToFit<Container, std::enable_if_t<
            std::is_same<decltype(std::declval<Container>().shrink_to_fit()), void>::value>>
    : std::true_type
    {};
    
    template <typename Container>
    void shrink_to_fit_if_enabled(Container& c, std::true_type)
    {
        c.shrink_to_fit();
    }
    
    template <typename Container>
    void shrink_to_fit_if_enabled(Container& c, std::false_type) {}
} // namespace detail

template <typename Container>
void reserve_if_enabled(Container& c, typename Container::size_type n)
{
    detail::reserve_if_enabled(c, n, detail::HasReserve<Container> {});
}

template <typename Container>
void shrink_to_fit_if_enabled(Container& c)
{
    detail::shrink_to_fit_if_enabled(c, detail::HasShrinkToFit<Container> {});
}

namespace detail
{
    template <typename T, typename = void>
    struct IsIterator
    {
        static constexpr bool value = false;
    };
    
    template <typename T>
    struct IsIterator<T, typename std::enable_if_t<
            !std::is_same<typename std::iterator_traits<T>::value_type, void>::value>>
    {
        static constexpr bool value = true;
    };
} // namespace detail

template <typename T>
constexpr bool is_iterator = detail::IsIterator<T>::value;

template <typename T, typename V = void>
using enable_if_iterator = std::enable_if_t<is_iterator<T>, V>;

namespace detail
{
    template <typename T, typename = void>
    struct IsMap : std::false_type {};
    
    template <typename T>
    struct IsMap<T, std::enable_if_t<
                        std::is_same<typename T::value_type,
                                std::pair<const typename T::key_type,
                                            typename T::mapped_type>
                        >::value>
                > : std::true_type {};
}

template <typename T>
constexpr bool is_map = detail::IsMap<T>::value;

namespace { template <bool> struct MapTagImpl {}; }
using MapTag    = MapTagImpl<true>;
using NonMapTag = MapTagImpl<false>;

template <typename T>
using MapTagType = MapTagImpl<is_map<T>>;

template <typename T, typename V = void>
using enable_if_map = std::enable_if_t<is_map<T>, V>;

template <typename T, typename V = void>
using enable_if_not_map = std::enable_if_t<!is_map<T>, V>;

template <typename Container, typename ConstIterator>
typename Container::iterator remove_constness(Container& c, ConstIterator it)
{
    return c.erase(it, it);
}

} // namespace octopus

#endif
