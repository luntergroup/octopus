// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef merge_transform_hpp
#define merge_transform_hpp

#include <algorithm>
#include <iterator>
#include <functional>

namespace octopus {

namespace detail {
    template <typename T1, typename T2 = T1>
    struct less
    {
        bool operator()(const T1& lhs, const T1& rhs) const {return lhs < rhs;}
        
        bool operator()(const T1& lhs, const T2& rhs) const {return lhs < rhs;}
        
        bool operator()(const T2& lhs, const T1& rhs) const {return lhs < rhs;}
        
        bool operator()(const T2& lhs, const T2& rhs) const {return lhs < rhs;}
    };
} // namespace detail

template <typename InputIt1, typename InputIt2, typename InputIt3, typename OutputIt,
typename BinaryOperation, typename UnaryOperation, typename Compare>
OutputIt merge_transform(InputIt1 first1, InputIt1 last1, InputIt2 first2,
                         InputIt3 first3, InputIt3 last3,
                         OutputIt d_first,
                         BinaryOperation binary_op, UnaryOperation unary_op,
                         Compare cmp)
{
    for (; first1 != last1; ++d_first) {
        if (first3 == last3) {
            return std::transform(first1, last1, first2, d_first, binary_op);
        }
        if (cmp(*first3, *first1)) {
            *d_first = unary_op(*first3);
            ++first3;
        } else {
            *d_first = binary_op(*first1, *first2);
            ++first1;
            ++first2;
        }
    }
    return std::transform(first3, last3, d_first, unary_op);
}

template <typename InputIt1, typename InputIt2, typename InputIt3, typename OutputIt,
typename BinaryOperation, typename UnaryOperation>
OutputIt merge_transform(InputIt1 first1, InputIt1 last1, InputIt2 first2,
                         InputIt3 first3, InputIt3 last3,
                         OutputIt d_first,
                         BinaryOperation binary_op, UnaryOperation unary_op)
{
    using V1 = typename std::iterator_traits<InputIt1>::value_type;
    using V2 = typename std::iterator_traits<InputIt3>::value_type;
    return merge_transform(first1, last1, first2, first3, last3, d_first, binary_op, unary_op,
                           detail::less<V1, V2>());
}

template <typename Container1, typename Container2, typename Container3, typename OutputIt,
typename BinaryOperation, typename UnaryOperation, typename Compare>
OutputIt merge_transform(const Container1& container1, const Container2& container2,
                         const Container3& container3,
                         OutputIt d_first,
                         BinaryOperation binary_op, UnaryOperation unary_op,
                         Compare cmp)
{
    return merge_transform(std::cbegin(container1), std::cend(container1), std::cbegin(container2),
                           std::cbegin(container3), std::cend(container3), d_first,
                           binary_op, unary_op, cmp);
}

template <typename Container1, typename Container2, typename Container3, typename OutputIt,
typename BinaryOperation, typename UnaryOperation>
OutputIt merge_transform(const Container1& container1, const Container2& container2,
                         const Container3& container3,
                         OutputIt d_first,
                         BinaryOperation binary_op, UnaryOperation unary_op)
{
    return merge_transform(std::cbegin(container1), std::cend(container1), std::cbegin(container2),
                           std::cbegin(container3), std::cend(container3), d_first,
                           binary_op, unary_op);
}

} // namespace octopus

#endif
