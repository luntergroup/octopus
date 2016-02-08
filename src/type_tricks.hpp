//
//  type_tricks.hpp
//  Octopus
//
//  Created by Daniel Cooke on 10/01/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#ifndef type_tricks_hpp
#define type_tricks_hpp

#include <type_traits>

namespace detail
{
    template <typename Container, typename = void>
    struct HasReserve : std::false_type {};
    
    template <typename Container>
    struct HasReserve<Container,
    std::enable_if_t<
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
    struct HasShrinkToFit<Container,
    std::enable_if_t<std::is_same<decltype(std::declval<Container>().shrink_to_fit()), void>::value>>
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
    struct IsIterator<T, typename std::enable_if_t<!std::is_same<typename std::iterator_traits<T>::value_type, void>::value>>
    {
        static constexpr bool value = true;
    };
} // namespace detail

template <typename T>
constexpr bool is_iterator = detail::IsIterator<T>::value;

template <typename T, typename V = void>
using enable_if_iterator = std::enable_if_t<is_iterator<T>, V>;

#endif /* type_tricks_hpp */
