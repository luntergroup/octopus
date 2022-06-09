// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef parallel_transform_hpp
#define parallel_transform_hpp

#include <iterator>
#include <vector>
#include <future>
#include <algorithm>
#include <cstddef>
#include <utility>
#include <type_traits>

#include <boost/optional.hpp>

#include "thread_pool.hpp"

namespace octopus {

namespace detail {

template <typename InputIt,
          typename OutputIt,
          typename UnaryOp>
OutputIt parallel_transform(InputIt first, InputIt last, OutputIt result, UnaryOp op,
                            std::random_access_iterator_tag)
{
    using value_type  = typename std::iterator_traits<InputIt>::value_type;
    using result_type = std::result_of_t<UnaryOp(value_type)>;
    std::vector<std::future<result_type>> results(std::distance(first, last));
    std::transform(first, last, std::begin(results),
                   [&op](const auto& value) {
                       return std::async(std::launch::async, op, std::cref(value));
                   });
    return std::transform(std::begin(results), std::end(results), result,
                          [](auto& f) { return f.get(); });
}

template <typename InputIt,
          typename OutputIt,
          typename UnaryOp>
OutputIt parallel_transform(InputIt first, InputIt last, OutputIt result, UnaryOp op,
                            std::input_iterator_tag)
{
    return std::transform(first, last, result, std::move(op));
}

template <typename InputIt1,
          typename InputIt2,
          typename OutputIt,
          typename BinaryOp>
OutputIt parallel_transform(InputIt1 first1, InputIt1 last1, InputIt2 first2, OutputIt result, BinaryOp op,
                            std::random_access_iterator_tag, std::random_access_iterator_tag)
{
    using value_type1  = typename std::iterator_traits<InputIt1>::value_type;
    using value_type2  = typename std::iterator_traits<InputIt2>::value_type;
    using result_type = std::result_of_t<BinaryOp(value_type1, value_type2)>;
    std::vector<std::future<result_type>> results(std::distance(first1, last1));
    std::transform(first1, last1, first2, std::begin(results),
                   [&op](const auto& a, const auto& b) {
                       return std::async(std::launch::async, op, std::cref(a), std::cref(b));
                   });
    return std::transform(std::begin(results), std::end(results), result,
                          [](auto& f) { return f.get(); });
}

template <typename InputIt1,
          typename InputIt2,
          typename OutputIt,
          typename BinaryOp>
OutputIt parallel_transform(InputIt1 first1, InputIt1 last1, InputIt2 first2, OutputIt result, BinaryOp op,
                            std::input_iterator_tag, std::input_iterator_tag)
{
    return std::transform(first1, last1, first2, result, std::move(op));
}

} // namespace detail

template <typename InputIt,
          typename OutputIt,
          typename UnaryOp>
OutputIt parallel_transform(InputIt first, InputIt last, OutputIt result, UnaryOp op)
{
    return detail::parallel_transform(first, last, result, std::move(op),
                                      typename std::iterator_traits<InputIt>::iterator_category {});
}

template <typename InputIt1,
          typename InputIt2,
          typename OutputIt,
          typename BinaryOp>
OutputIt parallel_transform(InputIt1 first1, InputIt1 last1, InputIt2 first2, OutputIt result, BinaryOp op)
{
    return detail::parallel_transform(first1, last1, first2, result, std::move(op),
                                      typename std::iterator_traits<InputIt1>::iterator_category {},
                                      typename std::iterator_traits<InputIt2>::iterator_category {});
}

namespace detail {

template <typename InputIt,
          typename OutputIt,
          typename UnaryOp>
OutputIt 
transform(InputIt first, InputIt last, OutputIt result, 
          UnaryOp op, ThreadPool& pool,
          std::random_access_iterator_tag)
{
    using value_type  = typename std::iterator_traits<InputIt>::value_type;
    using result_type = std::result_of_t<UnaryOp(value_type)>;
    const auto n = std::distance(first, last);
    if (n < 2) {
        return std::transform(first, last, result, std::move(op));
    }
    std::vector<std::future<result_type>> results(n);
    std::transform(first, std::prev(last), std::begin(results),
                   [&op, &pool] (auto&& value) {
                       return pool.try_push([&] () { return op(value); });
                   });
    // run last input in calling thread
    std::promise<result_type> last_result {};
    results.back() = last_result.get_future();
    last_result.set_value(op(*std::prev(last)));
    return std::transform(std::begin(results), std::end(results), result,
                          [](auto& f) { return f.get(); });
}

template <typename InputIt,
          typename OutputIt,
          typename UnaryOp>
OutputIt transform(InputIt first, InputIt last, OutputIt result, UnaryOp op, ThreadPool& pool,
                   std::input_iterator_tag)
{
    return std::transform(first, last, result, std::move(op));
}

template <typename InputIt1,
          typename InputIt2,
          typename OutputIt,
          typename BinaryOp>
OutputIt 
transform(InputIt1 first1, InputIt1 last1, InputIt2 first2, OutputIt result, 
          BinaryOp op, ThreadPool& pool,
          std::random_access_iterator_tag, std::random_access_iterator_tag)
{
    using value_type1 = typename std::iterator_traits<InputIt1>::value_type;
    using value_type2 = typename std::iterator_traits<InputIt2>::value_type;
    using result_type = std::result_of_t<BinaryOp(value_type1, value_type2)>;
    const auto n = std::distance(first1, last1);
    if (n < 2) {
        return std::transform(first1, last1, first2, result, std::move(op));
    }
    std::vector<std::future<result_type>> results(n);
    std::transform(first1, std::prev(last1), first2, std::begin(results),
                   [&op, &pool] (auto&& a, auto&& b) {
                       return pool.try_push([&] () { return op(a, b); });
                   });
    // run last input in calling thread
    std::promise<result_type> last_result {};
    results.back() = last_result.get_future();
    last_result.set_value(op(*std::prev(last1), *std::next(first2, n - 1)));
    return std::transform(std::begin(results), std::end(results), result,
                          [](auto& f) { return f.get(); });
}

template <typename InputIt1,
          typename InputIt2,
          typename OutputIt,
          typename BinaryOp>
OutputIt transform(InputIt1 first1, InputIt1 last1, InputIt2 first2, OutputIt result, BinaryOp op, ThreadPool& pool,
                   std::input_iterator_tag, std::input_iterator_tag)
{
    return std::transform(first1, last1, first2, result, std::move(op));
}

} // namespace detail

template <typename InputIt,
          typename OutputIt,
          typename UnaryOp>
OutputIt transform(InputIt first, InputIt last, OutputIt result, UnaryOp op, ThreadPool& pool)
{
    if (pool.size() > 1) {
        return detail::transform(first, last, result, std::move(op), pool,
                                 typename std::iterator_traits<InputIt>::iterator_category {});
    } else {
        return std::transform(first, last, result, std::move(op));
    }
}

template <typename InputIt,
          typename OutputIt,
          typename UnaryOp>
OutputIt transform(InputIt first, InputIt last, OutputIt result, UnaryOp op, boost::optional<ThreadPool&> pool)
{
    if (pool) {
        return transform(first, last, result, std::move(op), *pool);
    } else {
        return std::transform(first, last, result, std::move(op));
    }
}

template <typename InputIt1,
          typename InputIt2,
          typename OutputIt,
          typename BinaryOp>
OutputIt transform(InputIt1 first1, InputIt1 last1, InputIt2 first2, OutputIt result, BinaryOp op, ThreadPool& pool)
{
    if (pool.size() > 1) {
        return detail::transform(first1, last1, first2, result, std::move(op), pool,
                             typename std::iterator_traits<InputIt1>::iterator_category {},
                             typename std::iterator_traits<InputIt2>::iterator_category {});
    } else {
        return std::transform(first1, last1, first2, result, std::move(op));
    }
    
}

template <typename InputIt1,
          typename InputIt2,
          typename OutputIt,
          typename BinaryOp>
OutputIt transform(InputIt1 first1, InputIt1 last1, InputIt2 first2, OutputIt result, BinaryOp op, boost::optional<ThreadPool&> pool)
{
    if (pool) {
        return transform(first1, last1, first2, result, std::move(op), *pool);
    } else {
        return std::transform(first1, last1, first2, result, std::move(op));
    }
}

// for_each

namespace detail {

template <typename InputIt,
          typename UnaryOp>
void for_each(InputIt first, InputIt last, UnaryOp op, ThreadPool& pool, 
              std::input_iterator_tag)
{
    if (first == last) return;
    std::vector<std::future<void>> futures {};
    std::transform(first, last, std::back_inserter(futures),
                   [&op, &pool] (auto& value) {
                       return pool.try_push([&] () { return op(value); });
                   });
    for (auto& f : futures) f.get();
}

template <typename InputIt,
          typename UnaryOp>
void for_each(InputIt first, InputIt last, UnaryOp op, ThreadPool& pool, 
              std::bidirectional_iterator_tag)
{
    const auto n = std::distance(first, last);
    if (n < 2) {
        std::for_each(first, last, std::move(op));
    }
    std::vector<std::future<void>> futures(n - 1);
    std::transform(first, std::prev(last), std::begin(futures),
                   [&op, &pool] (auto& value) {
                       return pool.try_push([&] () { return op(value); });
                   });
    // run last input in calling thread
    op(*std::prev(last));
    for (auto& f : futures) f.get();
}

} // namespace detail

template <typename InputIt,
          typename UnaryOp>
void for_each(InputIt first, InputIt last, UnaryOp op, ThreadPool& pool)
{
    if (pool.size() > 1) {
        detail::for_each(first, last, std::move(op), pool);
    } else {
        std::for_each(first, last, std::move(op));
    }
}

template <typename InputIt,
          typename UnaryOp>
void for_each(InputIt first, InputIt last, UnaryOp op, boost::optional<ThreadPool&> pool)
{
    if (pool) {
        for_each(first, last, std::move(op), *pool);
    } else {
        std::for_each(first, last, std::move(op));
    }
}

} // namespace octopus

#endif
