//
//  parallel_transform.hpp
//
//  Created by Daniel Cooke on 29/10/2014.
//  Copyright (c) 2014 Oxford University. All rights reserved.
//

#ifndef parallel_transform_hpp
#define parallel_transform_hpp

#include <iterator>
#include <vector>
#include <future>
#include <algorithm>
#include <cstddef>
#include <utility>
#include <type_traits>

class sequential_execution_policy {};
extern const sequential_execution_policy seq;
class parallel_execution_policy {};
extern const parallel_execution_policy par;
class vector_execution_policy {};
extern const vector_execution_policy vec;

namespace detail
{
    template <typename ExecutionPolicy,
              typename InputIt,
              typename OutputIt,
              typename UnaryOp>
    OutputIt transform(ExecutionPolicy&& policy,
                       InputIt first, InputIt last, OutputIt result,
                       const UnaryOp unary_op,
                       std::random_access_iterator_tag)
    {
        using value_type  = typename std::iterator_traits<InputIt>::value_type;
        using result_type = std::result_of_t<UnaryOp(value_type)>;
        
        std::vector<std::future<result_type>> results(std::distance(first, last));
        
        std::transform(first, last, std::begin(results),
                       [&unary_op] (const auto& value) {
                           return std::async(unary_op, value);
                       });
        
        return std::transform(std::begin(results), std::end(results), result,
                              [] (auto& f) { return f.get(); });
    }
    
    template <typename ExecutionPolicy,
              typename InputIt,
              typename OutputIt,
              typename UnaryOp>
    OutputIt transform(ExecutionPolicy&& policy,
                       InputIt first, InputIt last, OutputIt result,
                       UnaryOp unary_op,
                       std::input_iterator_tag)
    {
        return std::transform(first, last, result, std::move(unary_op));
    }
    
    template <typename ExecutionPolicy,
              typename InputIt1,
              typename InputIt2,
              typename OutputIt,
              typename BinaryOp>
    OutputIt transform(ExecutionPolicy&& policy,
                       InputIt1 first1, InputIt1 last1, InputIt2 first2, OutputIt result,
                       BinaryOp unary_op,
                       std::random_access_iterator_tag, std::random_access_iterator_tag)
    {
        using value_type1  = typename std::iterator_traits<InputIt1>::value_type;
        using value_type2  = typename std::iterator_traits<InputIt2>::value_type;
        using result_type = std::result_of_t<BinaryOp(value_type1, value_type2)>;
        
        std::vector<std::future<result_type>> results(std::distance(first1, last1));
        
        std::transform(first1, last1, first2, std::begin(results),
                       [&unary_op] (const auto& a, const auto& b) {
                           return std::async(unary_op, a, b);
                       });
        
        return std::transform(std::begin(results), std::end(results), result,
                              [] (auto& f) { return f.get(); });
    }
    
    template <typename ExecutionPolicy,
              typename InputIt1,
              typename InputIt2,
              typename OutputIt,
              typename BinaryOp>
    OutputIt transform(ExecutionPolicy&& policy,
                       InputIt1 first1, InputIt1 last1, InputIt2 first2, OutputIt result,
                       BinaryOp unary_op,
                       std::input_iterator_tag, std::input_iterator_tag)
    {
        return std::transform(first1, last1, first2, result, std::move(unary_op));
    }
} // namespace detail

template <typename ExecutionPolicy,
          typename InputIt,
          typename OutputIt,
          typename UnaryOp>
OutputIt transform(ExecutionPolicy&& policy,
                   InputIt first, InputIt last, OutputIt result,
                   UnaryOp unary_op)
{
    return detail::transform(std::forward<ExecutionPolicy>(policy), first, last, result,
                             std::move(unary_op),
                             typename std::iterator_traits<InputIt>::iterator_category {});
}

template <typename ExecutionPolicy,
          typename InputIt1,
          typename InputIt2,
          typename OutputIt,
          typename BinaryOp>
OutputIt transform(ExecutionPolicy&& policy,
                   InputIt1 first1, InputIt1 last1, InputIt2 first2, OutputIt result,
                   BinaryOp unary_op)
{
    return detail::transform(std::forward<ExecutionPolicy>(policy), first1, last1, first2, result,
                             std::move(unary_op),
                             typename std::iterator_traits<InputIt1>::iterator_category {},
                             typename std::iterator_traits<InputIt2>::iterator_category {});
}

#endif
