// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef maths_hpp
#define maths_hpp

#include <vector>
#include <cstddef>
#include <cmath>
#include <numeric>
#include <algorithm>
#include <type_traits>
#include <iterator>
#include <functional>
#include <limits>
#include <utility>
#include <stdexcept>
#include <cassert>

#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/special_functions/sign.hpp>
#include <boost/math/special_functions/digamma.hpp>
#include <boost/math/special_functions/beta.hpp>
#include <boost/math/distributions/beta.hpp>
#include <boost/math/distributions/geometric.hpp>
#include <boost/math/distributions/binomial.hpp>
#include <boost/math/distributions/normal.hpp>

#include "htslib/kfunc.h"

#if defined(__SSE2__)
#include "fmath.hpp"
#endif // defined(__SSE2__)

namespace octopus { namespace maths {

namespace constants
{
    template <typename T = double>
    constexpr T ln10Div10 = T {0.230258509299404568401799145468436420760110148862877297603};
}

template <typename RealType>
bool is_subnormal(const RealType x) noexcept
{
    return std::fpclassify(x) == FP_SUBNORMAL;
}

template <typename RealType,
          typename = std::enable_if_t<std::is_floating_point<RealType>::value>>
RealType round(const RealType val, const unsigned precision = 2)
{
    const auto factor = std::pow(RealType {10.0}, precision);
    return std::round(val * factor) / factor;
}

template <typename T, typename = typename std::enable_if_t<std::is_floating_point<T>::value>>
T round_sf(const T x, const int n)
{
    // https://stackoverflow.com/a/13094362/2970186
    if (x == 0.0) return 0;
    auto factor = std::pow(10.0, n - std::ceil(std::log10(std::abs(x))));
    return std::round(x * factor) / factor;
}

template <typename T, typename = typename std::enable_if_t<std::is_floating_point<T>::value>>
bool almost_equal(const T lhs, T rhs, const int ulp = 1)
{
    return lhs == rhs || std::abs(lhs - rhs) < std::numeric_limits<T>::epsilon() * std::abs(lhs + rhs) * ulp;
}

template <typename T, typename = typename std::enable_if_t<std::is_floating_point<T>::value>>
bool almost_zero(const T x, const int ulp = 1)
{
    return almost_equal(x, T {0}, ulp);
}

template <typename T, typename = typename std::enable_if_t<std::is_floating_point<T>::value>>
bool almost_one(const T x, const int ulp = 1)
{
    return almost_equal(x, T {1}, ulp);
}

template <typename T, typename = typename std::enable_if_t<std::is_floating_point<T>::value>>
int count_leading_zeros(const T x)
{
    if (x == 0.0) return 0;
    return -std::ceil(std::log10(std::abs(x - std::numeric_limits<T>::epsilon())));
}

template <typename IntegerType,
          typename = std::enable_if_t<std::is_integral<IntegerType>::value>>
constexpr IntegerType ipow(const IntegerType base, const IntegerType exponent) noexcept
{
    if (exponent == 0) return 1;
    if (exponent == 1) return base;
    const auto y = ipow(base, exponent / 2);
    return exponent % 2 == 0 ? y * y : base * y * y;
}

template <typename IntegerType,
          typename = std::enable_if_t<std::is_integral<IntegerType>::value>>
bool is_safe_ipow(const IntegerType base, const IntegerType exponent) noexcept
{
    return exponent * std::log(base) <= std::log(std::numeric_limits<IntegerType>::max());
}

template <typename RealType>
constexpr RealType exp_maclaurin(const RealType x)
{
    return (6 + x * (6 + x * (3 + x))) * 0.16666666;
}

template <typename RealType>
constexpr RealType mercator(const RealType x)
{
    return x - x * x / 2 + x * x * x / 3;
}

struct IdFunction
{
    template <typename T>
    const T& operator()(const T& x) const noexcept { return x; }
};

template <typename RealType = double, typename InputIt, typename UnaryOperation>
auto mean(InputIt first, InputIt last, UnaryOperation unary_op)
{
    return std::accumulate(first, last, RealType {0},
                           [&] (const auto curr, const auto& x) {
                               return curr + unary_op(x);
                           }) / std::distance(first, last);
}

template <typename InputIt>
auto mean(InputIt first, InputIt last)
{
    return mean(first, last, IdFunction {});
}

template <typename Container>
auto mean(const Container& values)
{
    return mean(std::cbegin(values), std::cend(values));
}

template <typename Container, typename UnaryOperation>
auto mean(const Container& values, UnaryOperation unary_op)
{
    return mean(std::cbegin(values), std::cend(values), unary_op);
}

namespace detail {

template <typename T = double, typename ForwardIt>
T median_unsorted(ForwardIt first, ForwardIt last)
{
    const auto n = std::distance(first, last);
    assert(n > 0);
    if (n == 1) return *first;
    if (n == 2) return static_cast<T>(*first + *std::next(first)) / 2;
    const auto middle = std::next(first, n / 2);
    std::nth_element(first, middle, last);
    if (n % 2 == 1) {
        return *middle;
    } else {
        auto prev_middle_itr = std::max_element(first, middle);
        return static_cast<T>(*prev_middle_itr + *middle) / 2;
    }
}

template <typename T = double, typename ForwardIt>
T median_sorted(ForwardIt first, ForwardIt last)
{
    const auto n = std::distance(first, last);
    assert(n > 0);
    if (n == 1) return *first;
    const auto middle = std::next(first, n / 2);
    if (n % 2 == 1) {
        return *middle;
    } else {
        return static_cast<T>(*std::prev(middle) + *middle) / 2;
    }
}

template <typename T = double, typename ForwardIt>
T median_const(ForwardIt first, ForwardIt last)
{
    if (std::is_sorted(first, last)) {
        return median_sorted<T>(first, last);
    } else {
        std::vector<typename std::iterator_traits<ForwardIt>::value_type> tmp {first, last};
        return median_unsorted<T>(std::begin(tmp), std::end(tmp));
    }
}

} // namespace detail

template <typename T = double, typename ForwardIt>
T median(ForwardIt first, ForwardIt last)
{
    return detail::median_unsorted<T>(first, last);
}

template <typename T = double, typename Range>
T median(Range& values)
{
    return median<T>(std::begin(values), std::end(values));
}

template <typename T = double, typename Range>
T median(const Range& values)
{
    return detail::median_const<T>(std::cbegin(values), std::cend(values));
}

template <typename ForwardIterator, typename UnaryOperation>
auto stdev(ForwardIterator first, ForwardIterator last, UnaryOperation unary_op)
{
    const auto m = mean(first, last, unary_op);
    const auto n = std::distance(first, last);
    const auto sum_square = [&] (auto total, const auto& x) { return total + std::pow(unary_op(x) - m, 2); };
    const auto ss = std::accumulate(first, last, 0.0, sum_square);
    return std::sqrt(ss / n);
}

template <typename ForwardIterator>
auto stdev(ForwardIterator first, ForwardIterator last)
{
    return stdev(first, last, IdFunction {});
}

template <typename Container>
auto stdev(const Container& values)
{
    return stdev(std::cbegin(values), std::cend(values));
}

template <typename Container, typename UnaryOperation>
auto stdev(const Container& values, UnaryOperation unary_op)
{
    return stdev(std::cbegin(values), std::cend(values), unary_op);
}

template <typename RealType = double, typename InputIt>
RealType rmq(InputIt first, InputIt last)
{
    if (first == last) return 0.0;
    return std::sqrt((std::inner_product(first, last, first, RealType {0}))
                     / static_cast<RealType>(std::distance(first, last)));
}

template <typename RealType = double, typename Container>
RealType rmq(const Container& values)
{
    return rmq<RealType>(std::cbegin(values), std::cend(values));
}

inline float fast_exp(const float x) noexcept
{
    #if defined(__SSE2__)
    return fmath::exp(x);
    #else
    return std::exp(x);
    #endif // defined(__SSE2__)
}

template <typename RealType,
          typename = std::enable_if_t<std::is_floating_point<RealType>::value>>
inline RealType fast_exp(const RealType x) noexcept
{
    return std::exp(x);
}

inline float fast_log(const float x) noexcept
{
    #if defined(__SSE2__)
    return fmath::log(x);
    #else
    return std::log(x);
    #endif // defined(__SSE2__)
}

template <typename RealType,
          typename = std::enable_if_t<std::is_floating_point<RealType>::value>>
inline RealType fast_log(const RealType x) noexcept
{
    return std::log(x);
}

template <typename RealType,
          typename = std::enable_if_t<std::is_floating_point<RealType>::value>>
inline RealType log_sum_exp(const RealType a, const RealType b)
{
    const auto r = std::minmax(a, b);
    return r.second + std::log1p(std::exp(r.first - r.second));
}

template <typename RealType,
          typename = std::enable_if_t<std::is_floating_point<RealType>::value>>
inline RealType log_sum_exp(const RealType a, const RealType b, const RealType c)
{
    const auto max = std::max({a, b, c});
    return max + std::log(std::exp(a - max) + std::exp(b - max) + std::exp(c - max));
}

template <typename RealType,
          typename = std::enable_if_t<std::is_floating_point<RealType>::value>>
inline RealType log_sum_exp(std::initializer_list<RealType> il)
{
    const auto max = std::max(il);
    return max + std::log(std::accumulate(std::cbegin(il), std::cend(il), RealType {0},
                                          [max] (const auto curr, const auto x) {
                                              return curr + std::exp(x - max);
                                          }));
}

template <typename ForwardIt,
          typename = std::enable_if_t<!std::is_floating_point<ForwardIt>::value>>
inline auto log_sum_exp(ForwardIt first, ForwardIt last)
{
    assert(first != last);
    using RealType = typename std::iterator_traits<ForwardIt>::value_type;
    const auto max = *std::max_element(first, last);
    return max + std::log(std::accumulate(first, last, RealType {0},
                                          [max] (const auto curr, const auto x) {
                                              return curr + std::exp(x - max);
                                          }));
}

template <typename RealType,
          typename = std::enable_if_t<std::is_floating_point<RealType>::value>>
inline RealType log_sum_exp(const std::array<RealType, 1>& logs)
{
    return logs[0];
}

template <typename RealType,
          typename = std::enable_if_t<std::is_floating_point<RealType>::value>>
inline RealType log_sum_exp(const std::array<RealType, 2>& logs)
{
    return log_sum_exp(logs[0], logs[1]);
}

template <typename RealType,
          typename = std::enable_if_t<std::is_floating_point<RealType>::value>>
inline RealType log_sum_exp(const std::array<RealType, 3>& logs)
{
    return log_sum_exp(logs[0], logs[1], logs[2]);
}

template <typename Container>
inline auto log_sum_exp(const Container& values)
{
    return log_sum_exp(std::cbegin(values), std::cend(values));
}

template <typename RealType,
          typename = std::enable_if_t<std::is_floating_point<RealType>::value>>
inline RealType fast_log_sum_exp(const RealType a, const RealType b) noexcept
{
    const auto r = std::minmax(a, b);
    return r.second + fast_log(RealType {1} + fast_exp(r.first - r.second));
}

template <typename RealType,
          typename = std::enable_if_t<std::is_floating_point<RealType>::value>>
inline RealType fast_log_sum_exp(const RealType a, const RealType b, const RealType c) noexcept
{
    const auto max = std::max({a, b, c});
    return max + fast_log(fast_exp(a - max) + fast_exp(b - max) + fast_exp(c - max));
}

template <typename RealType,
          typename = std::enable_if_t<std::is_floating_point<RealType>::value>>
inline RealType fast_log_sum_exp(std::initializer_list<RealType> il) noexcept
{
    const auto max = std::max(il);
    return max + fast_log(std::accumulate(std::cbegin(il), std::cend(il), RealType {0},
                                          [max] (const auto curr, const auto x) noexcept {
                                              return curr + fast_exp(x - max);
                                          }));
}

template <typename ForwardIt,
          typename = std::enable_if_t<!std::is_floating_point<ForwardIt>::value>>
inline auto fast_log_sum_exp(ForwardIt first, ForwardIt last) noexcept
{
    assert(first != last);
    using RealType = typename std::iterator_traits<ForwardIt>::value_type;
    const auto max = *std::max_element(first, last);
    return max + fast_log(std::accumulate(first, last, RealType {0},
                                          [max] (const auto curr, const auto x) noexcept {
                                              return curr + fast_exp(x - max);
                                          }));
}

template <typename RealType,
          typename = std::enable_if_t<std::is_floating_point<RealType>::value>>
inline RealType fast_log_sum_exp(const std::array<RealType, 1>& logs) noexcept
{
    return logs[0];
}

template <typename RealType,
          typename = std::enable_if_t<std::is_floating_point<RealType>::value>>
inline RealType fast_log_sum_exp(const std::array<RealType, 2>& logs) noexcept
{
    return fast_log_sum_exp(logs[0], logs[1]);
}

template <typename RealType,
          typename = std::enable_if_t<std::is_floating_point<RealType>::value>>
inline RealType fast_log_sum_exp(const std::array<RealType, 3>& logs) noexcept
{
    return fast_log_sum_exp(logs[0], logs[1], logs[2]);
}

template <typename Container>
inline auto fast_log_sum_exp(const Container& values) noexcept
{
    return fast_log_sum_exp(std::cbegin(values), std::cend(values));
}

template <typename T, typename IntegerType,
          typename = std::enable_if_t<std::is_integral<IntegerType>::value>>
T factorial(const IntegerType x)
{
    return boost::math::factorial<double>(x);
}

template <typename RealType, typename IntegerType,
          typename = std::enable_if_t<std::is_floating_point<RealType>::value>,
          typename = std::enable_if_t<std::is_integral<IntegerType>::value>>
RealType log_factorial(IntegerType x)
{
    return std::lgamma(x + 1);
}

template <typename InputIt, typename Log>
auto entropy(InputIt first, InputIt last, Log&& log)
{
    using RealType = typename std::iterator_traits<InputIt>::value_type;
    static_assert(std::is_floating_point<RealType>::value,
                  "entropy is only defined for floating point values");
    const auto add_entropy = [&] (auto total, auto p) { return total + (p > 0 ? p * log(p) : 0); };
    return -std::accumulate(first, last, RealType {0}, add_entropy);
}

template <typename Range, typename Log>
auto entropy(const Range& values, Log&& log)
{
    return entropy(std::cbegin(values), std::cend(values), std::forward<Log>(log));
}

template <typename InputIt>
auto entropy(InputIt first, InputIt last)
{
    return entropy(first, last, [] (auto x) { return std::log(x); });
}

template <typename Range>
auto entropy(const Range& values)
{
    return entropy(std::cbegin(values), std::cend(values));
}

template <typename InputIt>
auto entropy2(InputIt first, InputIt last)
{
    return entropy(first, last, [] (auto x) { return std::log2(x); });
}

template <typename Range>
auto entropy2(const Range& values)
{
    return entropy2(std::cbegin(values), std::cend(values));
}

template <typename InputIt>
auto entropy10(InputIt first, InputIt last)
{
    return entropy(first, last, [] (auto x) { return std::log10(x); });
}

template <typename Range>
auto entropy10(const Range& values)
{
    return entropy10(std::cbegin(values), std::cend(values));
}

template <typename RealType, typename IntegerType,
          typename = std::enable_if_t<std::is_integral<IntegerType>::value>,
          typename = std::enable_if_t<std::is_floating_point<RealType>::value>>
RealType binomial_coefficient(const IntegerType n, const IntegerType k)
{
    return boost::math::binomial_coefficient<RealType>(n, k);
}

template <typename RealType, typename IntegerType,
          typename = std::enable_if_t<std::is_integral<IntegerType>::value>,
          typename = std::enable_if_t<std::is_floating_point<RealType>::value>>
RealType log_binomial_coefficient(const IntegerType n, const IntegerType k)
{
    return log_factorial<RealType>(n) - (log_factorial<RealType>(k) + log_factorial<RealType>(n - k));
}

template <typename IntegerType, typename RealType = double,
          typename = std::enable_if_t<std::is_integral<IntegerType>::value>,
          typename = std::enable_if_t<std::is_floating_point<RealType>::value>>
RealType fisher_exact_test(const IntegerType a, const IntegerType b, const IntegerType c, const IntegerType d)
{
    double fisher_left_p, fisher_right_p, fisher_twosided_p;
    kt_fisher_exact(a, b, c, d, &fisher_left_p, &fisher_right_p, &fisher_twosided_p);
    return fisher_twosided_p;
}

template <typename IntegerType, typename RealType,
          typename = std::enable_if_t<std::is_integral<IntegerType>::value>,
          typename = std::enable_if_t<std::is_floating_point<RealType>::value>>
RealType geometric_pdf(const IntegerType k, const RealType p)
{
    boost::math::geometric_distribution<RealType> dist {p};
    return boost::math::pdf(dist, k);
}

template <typename IntegerType, typename RealType,
          typename = std::enable_if_t<std::is_integral<IntegerType>::value>,
          typename = std::enable_if_t<std::is_floating_point<RealType>::value>>
RealType binomial_pdf(const IntegerType k, const IntegerType n, const RealType p)
{
    boost::math::binomial_distribution<RealType> dist {static_cast<RealType>(n), p};
    return boost::math::pdf(dist, k);
}

template <typename IntegerType, typename RealType,
          typename = std::enable_if_t<std::is_integral<IntegerType>::value>,
          typename = std::enable_if_t<std::is_floating_point<RealType>::value>>
RealType log_poisson_pmf(const IntegerType k, const RealType mu)
{
    if (k > 0) {
        return k * std::log(mu) - std::lgamma(k) - mu;
    } else {
        return -mu;
    }
}

template <typename IntegerType, typename RealType,
          typename = std::enable_if_t<std::is_integral<IntegerType>::value>,
          typename = std::enable_if_t<std::is_floating_point<RealType>::value>>
RealType poisson_cdf(const IntegerType k, const RealType mu)
{
    return almost_zero(mu) ? 1.0 : boost::math::gamma_q(k + 1, mu);
}

template <typename IntegerType, typename RealType,
          typename = std::enable_if_t<std::is_integral<IntegerType>::value>,
          typename = std::enable_if_t<std::is_floating_point<RealType>::value>>
RealType poisson_sf(const IntegerType k, const RealType mu)
{
    return almost_zero(mu) ? 0.0 : boost::math::gamma_p(k + 1, mu);
}

template <typename IntegerType, typename RealType,
          typename = std::enable_if_t<std::is_integral<IntegerType>::value>,
          typename = std::enable_if_t<std::is_floating_point<RealType>::value>>
RealType log_poisson_cdf(const IntegerType k, const RealType mu)
{
    return std::log(poisson_cdf(k, mu));
}

template <typename IntegerType, typename RealType,
          typename = std::enable_if_t<std::is_integral<IntegerType>::value>,
          typename = std::enable_if_t<std::is_floating_point<RealType>::value>>
RealType log_poisson_sf(const IntegerType k, const RealType mu)
{
    return std::log(poisson_sf(k, mu));
}

template <typename RealType,
          typename = std::enable_if_t<std::is_floating_point<RealType>::value>>
RealType normal_cdf(const RealType x, const RealType mu, const RealType sigma)
{
    boost::math::normal_distribution<RealType> dist {mu, sigma};
    return boost::math::cdf(dist, x);
}

template <typename RealType,
          typename = std::enable_if_t<std::is_floating_point<RealType>::value>>
RealType normal_sf(const RealType x, const RealType mu, const RealType sigma)
{
    boost::math::normal_distribution<RealType> dist {mu, sigma};
    return boost::math::cdf(boost::math::complement(dist, x));
}

template <typename ForwardIt>
auto log_beta(const ForwardIt first, const ForwardIt last)
{
    using T = std::decay_t<typename std::iterator_traits<ForwardIt>::value_type>;
    static_assert(std::is_floating_point<T>::value,
                  "log_beta is only defined for floating point types.");
    return std::accumulate(first, last, T {0}, [] (const auto curr, const auto x) { return curr + std::lgamma(x); })
           - std::lgamma(std::accumulate(first, last, T {0}));
}

template <typename Container>
auto log_beta(const Container& values)
{
    return log_beta(std::cbegin(values), std::cend(values));
}

template <typename ForwardIt1, typename ForwardIt2>
auto log_dirichlet(ForwardIt1 firstalpha, ForwardIt1 lastalpha, ForwardIt2 firstpi)
{
    using T = std::decay_t<typename std::iterator_traits<ForwardIt1>::value_type>;
    static_assert(std::is_floating_point<T>::value,
                  "log_dirichlet is only defined for floating point types.");
    return std::inner_product(firstalpha, lastalpha, firstpi, T {0}, std::plus<> {},
                              [] (const auto a, const auto p) { return (a - 1) * std::log(p); })
            - log_beta(firstalpha, lastalpha);
}

template <typename Container1, typename Container2>
auto log_dirichlet(const Container1& alpha, const Container2& pi)
{
    return log_dirichlet(std::cbegin(alpha), std::cend(alpha), std::cbegin(pi));
}

template <typename ForwardIt>
auto dirichlet_expectation(ForwardIt first_alpha, ForwardIt last_alpha)
{
    using T = std::decay_t<typename std::iterator_traits<ForwardIt>::value_type>;
    static_assert(std::is_floating_point<T>::value,
                  "log_dirichlet is only defined for floating point types.");
    const auto K = static_cast<std::size_t>(std::distance(first_alpha, last_alpha));
    const auto a0 = std::accumulate(first_alpha, last_alpha, T {0});
    std::vector<T> result(K);
    std::transform(first_alpha, last_alpha, std::begin(result), [a0] (auto a) { return a / a0; });
    return result;
}

template <typename Range>
auto dirichlet_expectation(const Range& values)
{
    return dirichlet_expectation(std::cbegin(values), std::cend(values));
}

template <typename ForwardIt>
auto dirichlet_expectation(const unsigned i, ForwardIt first_alpha, ForwardIt last_alpha)
{
    using T = std::decay_t<typename std::iterator_traits<ForwardIt>::value_type>;
    static_assert(std::is_floating_point<T>::value,
                  "log_dirichlet is only defined for floating point types.");
    assert(i < static_cast<unsigned>(std::distance(first_alpha, last_alpha)));
    return *std::next(first_alpha, i) / std::accumulate(first_alpha, last_alpha, T {0});
}

template <typename Range>
auto dirichlet_expectation(const unsigned i, const Range& values)
{
    return dirichlet_expectation(i, std::cbegin(values), std::cend(values));
}

template <typename ForwardIt>
auto dirichlet_entropy(ForwardIt first_alpha, ForwardIt last_alpha)
{
    using T = std::decay_t<typename std::iterator_traits<ForwardIt>::value_type>;
    static_assert(std::is_floating_point<T>::value,
                  "log_dirichlet is only defined for floating point types.");
    const auto K = static_cast<T>(std::distance(first_alpha, last_alpha));
    const auto a0 = std::accumulate(first_alpha, last_alpha, T {0});
    using boost::math::digamma;
    return log_beta(first_alpha, last_alpha) + (a0 - K) * digamma(a0)
           - std::accumulate(first_alpha, last_alpha, T {0}, [] (auto curr, auto a) { return curr + (a - 1) * digamma(a); });
}

template <typename Range>
auto dirichlet_entropy(const Range& values)
{
    return dirichlet_entropy(std::cbegin(values), std::cend(values));
}

template <typename RealType, typename Iterator>
inline RealType log_multinomial_coefficient(Iterator first, Iterator last)
{
    using IntegerType = typename std::iterator_traits<Iterator>::value_type;
    return log_factorial<RealType>(std::accumulate(first, last, IntegerType {0}))
        - std::accumulate(first, last, RealType {0}, [] (auto total, IntegerType x) { 
            return total + log_factorial<RealType>(x); });
}

template <typename RealType, typename IntegerType>
inline RealType log_multinomial_coefficient(std::initializer_list<IntegerType> il)
{
    return log_multinomial_coefficient(std::cbegin(il), std::cend(il));
}

template <typename RealType, typename Container>
inline RealType log_multinomial_coefficient(const Container& values)
{
    return log_multinomial_coefficient<RealType>(std::cbegin(values), std::cend(values));
}

template <typename RealType, typename IntegerType>
inline IntegerType multinomial_coefficient(std::initializer_list<IntegerType> il)
{
    return static_cast<IntegerType>(std::exp(log_multinomial_coefficient<RealType, IntegerType>(std::move(il))));
}

template <typename IntegerType, typename RealType, typename Iterator>
inline IntegerType multinomial_coefficient(Iterator first, Iterator last)
{
    return static_cast<IntegerType>(std::exp(log_multinomial_coefficient<RealType>(first, last)));
}

template <typename IntegerType, typename RealType, typename Container>
inline IntegerType multinomial_coefficient(const Container& values)
{
    return multinomial_coefficient<IntegerType, RealType>(std::cbegin(values), std::cend(values));
}

template <typename RealType, typename ForwardIt1, typename ForwardIt2>
inline RealType multinomial_pdf(ForwardIt1 first_z, ForwardIt1 last_z, ForwardIt2 first_p)
{
    auto r = std::inner_product(first_z, last_z, first_p, RealType {0}, std::multiplies<> {},
                                [] (auto z_i, auto p_i) { return std::pow(p_i, z_i); });
    return multinomial_coefficient<RealType>(first_z, last_z) * r;
}

template <typename IntegerType, typename RealType>
inline RealType multinomial_pdf(const std::vector<IntegerType>& z, const std::vector<RealType>& p)
{
    assert(z.size() == p.size());
    return multinomial_pdf<RealType>(std::cbegin(z), std::cend(z), std::cbegin(p));
}

template <typename RealType, typename ForwardIt1, typename ForwardIt2>
inline RealType log_multinomial_pdf(ForwardIt1 first_z, ForwardIt1 last_z, ForwardIt2 first_p)
{
    auto r = std::inner_product(first_z, last_z, first_p, RealType {0}, std::plus<> {},
                                [] (auto z_i, auto p_i) { return z_i > 0 ? z_i * std::log(p_i) : 0.0; });
    return log_multinomial_coefficient<RealType>(first_z, last_z) + r;
}

template <typename IntegerType, typename RealType>
inline RealType log_multinomial_pdf(const std::vector<IntegerType>& z, const std::vector<RealType>& p)
{
    assert(z.size() == p.size());
    return log_multinomial_pdf<RealType>(std::cbegin(z), std::cend(z), std::cbegin(p));
}

// Returns approximate y such that digamma(y) = x
template <typename RealType>
inline RealType digamma_inv(const RealType x, const RealType epsilon = 10e-8)
{
    RealType l {1};
    auto y = std::exp(x);
    while (l > epsilon) {
        y += l * boost::math::sign(x - boost::math::digamma<RealType>(y));
        l /= 2;
    }
    return y;
}

namespace detail {

template <typename T, typename RealType>
T ifactorial(RealType x, std::true_type)
{
    return factorial<T, unsigned>(x);
}

template <typename T, typename RealType>
T ifactorial(RealType x, std::false_type)
{
    return factorial<T>(x);
}

template <typename T, typename RealType>
T ifactorial(RealType x)
{
    return ifactorial<T>(x, std::is_floating_point<RealType> {});
}

} // namespace detail

template <typename RealType>
RealType dirichlet_multinomial(const RealType z1, const RealType z2, const RealType a1, const RealType a2)
{
    auto z_0 = z1 + z2;
    auto a_0 = a1 + a2;
    using detail::ifactorial;
    auto z_m = ifactorial<RealType>(z1) * ifactorial<RealType>(z2);
    return (ifactorial<RealType>(z_0) / z_m) *
            (std::tgamma(a_0) / std::tgamma(z_0 + a_0)) *
            (std::tgamma(z1 + a1) * std::tgamma(z2 + a2)) / (std::tgamma(a1) + std::tgamma(a2));
}

template <typename RealType>
RealType dirichlet_multinomial(const RealType z1, const RealType z2, const RealType z3,
                               const RealType a1, const RealType a2, const RealType a3)
{
    auto z_0 = z1 + z2 + z3;
    auto a_0 = a1 + a2 + a3;
    using detail::ifactorial;
    auto z_m = ifactorial<RealType>(z1) * ifactorial<RealType>(z2) * ifactorial<RealType>(z3);
    return (ifactorial<RealType>(z_0) / z_m) *
            (std::tgamma(a_0) / std::tgamma(z_0 + a_0)) *
            (std::tgamma(z1 + a1) * std::tgamma(z2 + a2) *
            std::tgamma(z3 + a3)) / (std::tgamma(a1) + std::tgamma(a2) + std::tgamma(a3));
}

template <typename RealType>
RealType dirichlet_multinomial(const std::vector<RealType>& z, const std::vector<RealType>& a)
{
    auto z_0 = std::accumulate(std::cbegin(z), std::cend(z), RealType {0});
    auto a_0 = std::accumulate(std::cbegin(a), std::cend(a), RealType {0});
    RealType z_m {1};
    using detail::ifactorial;
    for (auto z_i : z) {
        z_m *= ifactorial<RealType>(z_i);
    }
    RealType g {1};
    for (std::size_t i {0}; i < z.size(); ++i) {
        g *= std::tgamma(z[i] + a[i]) / std::tgamma(a[i]);
    }
    return (ifactorial<RealType>(z_0) / z_m) * (std::tgamma(a_0) / std::tgamma(z_0 + a_0)) * g;
}

template <typename RealType>
RealType beta_binomial(const RealType k, const RealType n, const RealType alpha, const RealType beta)
{
    return dirichlet_multinomial<RealType>(k, n - k, alpha, beta);
}

namespace detail {

template <typename RealType>
bool is_mldp_converged(std::vector<RealType>& lhs, const std::vector<RealType>& rhs,
                       const RealType epsilon)
{
    std::transform(std::cbegin(lhs), std::cend(lhs), std::cbegin(rhs), std::begin(lhs),
                   [] (const auto a, const auto b) { return std::abs(a - b); });
    return std::all_of(std::cbegin(lhs), std::cend(lhs),
                       [epsilon] (const auto x) { return x < epsilon; });
}

} // namespace detail

template <typename RealType>
std::vector<RealType>
dirichlet_mle(std::vector<RealType> pi, const RealType precision,
              const unsigned max_iterations = 100, const RealType epsilon = 0.0001)
{
    std::transform(std::cbegin(pi), std::cend(pi), std::begin(pi),
                   [] (const auto p) { return std::log(p); });
    const auto l = pi.size();
    const RealType u {RealType {1} / l};
    std::vector<RealType> result(l, u), curr_result(l, u), means(l, u);
    for (unsigned n {0}; n < max_iterations; ++n) {
        RealType v {0};
        for (std::size_t j {0}; j < l; ++j) {
            v += means[j] * (pi[j] - boost::math::digamma<RealType>(precision * means[j]));
        }
        for (std::size_t k {0}; k < l; ++k) {
            curr_result[k] = digamma_inv<RealType>(pi[k] - v);
            means[k] = curr_result[k] / std::accumulate(std::cbegin(curr_result), std::cend(curr_result), RealType {0});
        }
        if (detail::is_mldp_converged(result, curr_result, epsilon)) {
            return curr_result;
        }
        result = curr_result;
    }
    return result;
}

template <typename NumericType = float, typename RealType,
          typename = std::enable_if_t<std::is_floating_point<RealType>::value>>
NumericType probability_true_to_phred(const RealType p)
{
    return NumericType{-10} * std::log10(std::max(RealType {1} - p, std::numeric_limits<RealType>::epsilon()));
}

template <typename NumericType = float, typename RealType,
          typename = std::enable_if_t<std::is_floating_point<RealType>::value>>
NumericType probability_true_to_phred(const RealType p, const unsigned precision)
{
    return round(static_cast<NumericType>(RealType {-10} * std::log10(std::max(RealType {1} - p, std::numeric_limits<RealType>::epsilon()))), precision);
}

template <typename RealType = double, typename NumericType>
RealType phred_to_probability(const NumericType phred)
{
    return RealType {1} - std::pow(RealType {10}, RealType {-1} * static_cast<RealType>(phred) / RealType {10});
}

template <typename MapType>
typename MapType::key_type sum_keys(const MapType& map)
{
    return std::accumulate(std::cbegin(map), std::cend(map), typename MapType::key_type {},
                           [] (const auto previous, const auto& p) { return previous + p.first; });
}

template <typename ResultType, typename MapType, typename UnaryOperation>
ResultType sum_keys(const MapType& map, UnaryOperation op)
{
    return std::accumulate(std::cbegin(map), std::cend(map), ResultType {0},
                           [op] (const auto previous, const auto& p) { return previous + op(p.first); });
}

template <typename MapType>
typename MapType::mapped_type sum_values(const MapType& map)
{
    return std::accumulate(std::cbegin(map), std::cend(map), typename MapType::mapped_type {},
                           [] (const auto previous, const auto& p) { return previous + p.second; });
}

template <typename ResultType, typename MapType, typename UnaryOperation>
ResultType sum_values(const MapType& map, UnaryOperation op)
{
    return std::accumulate(std::cbegin(map), std::cend(map), ResultType {0},
                           [op] (const auto previous, const auto& p) { return previous + op(p.second); });
}

template <typename Map>
std::size_t sum_sizes(const Map& map)
{
    return std::accumulate(std::cbegin(map), std::cend(map), std::size_t {0},
                           [] (const auto& p, const auto& v) { return p + v.second.size(); });
}

template <typename InputIt1, typename InputIt2, typename InputIt3, typename T,
          typename BinaryOperation1, typename BinaryOperation2>
T inner_product(InputIt1 first1, InputIt1 last1,
                InputIt2 first2, InputIt3 first3, T value,
                BinaryOperation1 op1, BinaryOperation2 op2)
{
    while (first1 != last1) {
        value = op1(value, op2(*first1, *first2, *first3));
        ++first1;
        ++first2;
        ++first3;
    }
    return value;
}

template <typename InputIt1, typename InputIt2, typename InputIt3, typename InputIt4,
          typename T, typename BinaryOperation1, typename BinaryOperation2>
T inner_product(InputIt1 first1, InputIt1 last1,
                InputIt2 first2, InputIt3 first3,
                InputIt4 first4, T value,
                BinaryOperation1 op1, BinaryOperation2 op2)
{
    while (first1 != last1) {
        value = op1(value, op2(*first1, *first2, *first3, *first4));
        ++first1;
        ++first2;
        ++first3;
        ++first4;
    }
    return value;
}

template <typename RealType>
RealType beta_cdf(const RealType a, const RealType b, const RealType x)
{
    const boost::math::beta_distribution<> beta_dist {a, b};
    return boost::math::cdf(beta_dist, x);
}

template <typename RealType>
RealType beta_sf(const RealType a, const RealType b, const RealType x)
{
    const boost::math::beta_distribution<> beta_dist {a, b};
    return boost::math::cdf(boost::math::complement(beta_dist, x));
}

template <typename RealType>
RealType beta_tail_probability(const RealType a, const RealType b, const RealType x)
{
    return beta_cdf(a, b, x) + beta_sf(a, b, RealType {1} - x);
}

namespace detail {

template <typename RealType>
std::pair<RealType, RealType>
uniform_hdi(const RealType mass)
{
    const auto x = RealType {0.5} - mass / 2;
    return std::make_pair(x, x + mass);
}

template <typename RealType>
std::pair<RealType, RealType>
beta_hdi_symmetric(const RealType a, const RealType mass)
{
    const auto x = boost::math::ibeta_inv(a, a, (RealType {1} - mass) / 2);
    return std::make_pair(x, RealType {1} - x);
}

template <typename RealType>
std::pair<RealType, RealType>
beta_hdi_unbounded_rhs(const RealType a, const RealType mass)
{
    // Reverse J shaped
    return std::make_pair(boost::math::ibeta_inv(a, RealType {1}, RealType {1} - mass), RealType {1});
}

template <typename RealType>
std::pair<RealType, RealType>
beta_hdi_unbounded_lhs(const RealType b, const RealType mass)
{
    // J shaped
    return std::make_pair(RealType {0}, boost::math::ibeta_inv(RealType {1}, b, mass));
}

template <typename RealType>
std::pair<RealType, RealType>
beta_hdi_skewed(const RealType a, const RealType b, const RealType mass)
{
    const auto c = (RealType {1} - mass) / 2;
    return std::make_pair(boost::math::ibeta_inv(a, b, c), boost::math::ibeta_inv(a, b, c + mass));
}

} // namespace detail

template <typename RealType>
std::pair<RealType, RealType>
beta_hdi(RealType a, RealType b, const RealType mass)
{
    static_assert(std::is_floating_point<RealType>::value, "beta_hdi only works for floating point types");
    if (mass < RealType {0} || mass > RealType {1}) {
        throw std::domain_error {"beta_hdi: given mass not in range [0, 1]"};
    }
    if (a <= RealType {0} || b <= RealType {0}) {
        throw std::domain_error {"beta_hdi: given non-positive parameter"};
    }
    if (mass == RealType {0}) {
        const auto mean = a / (a + b);
        return std::make_pair(mean, mean);
    }
    if (mass == RealType {1}) {
        return std::make_pair(RealType {0}, RealType {1});
    }
    if (a == b) {
        if (a == RealType {1}) {
            return detail::uniform_hdi(mass);
        } else {
            return detail::beta_hdi_symmetric(a, mass);
        }
    }
    if (a == RealType {1}) {
        return detail::beta_hdi_unbounded_lhs(b, mass);
    }
    if (b == RealType {1}) {
        return detail::beta_hdi_unbounded_rhs(a, mass);
    }
    return detail::beta_hdi_skewed(a, b, mass);
}

template <typename RealType>
RealType dirichlet_variance(const std::vector<RealType>& alphas, const std::size_t k)
{
    const auto a_0 = std::accumulate(std::cbegin(alphas), std::cend(alphas), RealType {});
    return (alphas[k] * (a_0 - alphas[k])) / (a_0 * a_0 * (a_0 + 1));
}

template <typename RealType>
RealType dirichlet_marginal_cdf(const std::vector<RealType>& alphas, const std::size_t k, const RealType x)
{
    const auto a_0 = std::accumulate(std::cbegin(alphas), std::cend(alphas), RealType {});
    return beta_cdf(alphas[k], a_0 - alphas[k], x);
}

template <typename RealType>
RealType dirichlet_marginal_sf(const std::vector<RealType>& alphas, const std::size_t k, const RealType x)
{
    const auto a_0 = std::accumulate(std::cbegin(alphas), std::cend(alphas), RealType {});
    return beta_sf(alphas[k], a_0 - alphas[k], x);
}

template <typename Range>
void log_each(Range& values)
{
    for (auto& v : values) v = std::log(v);
}

template <typename Container>
Container log_each_copy(Container values)
{
    log_each(values);
    return values;
}

template <typename Range>
void exp_each(Range& values)
{
    for (auto& v : values) v = std::exp(v);
}

template <typename Container>
Container exp_each_copy(Container values)
{
    exp_each(values);
    return values;
}

template <typename ForwardIterator>
auto normalise(const ForwardIterator first, const ForwardIterator last)
{
    using T = typename std::iterator_traits<ForwardIterator>::value_type;
    const auto norm = std::accumulate(first, last, T {});
    if (norm > 0) std::for_each(first, last, [norm] (auto& value) { value /= norm; });
    return norm;
}

template <typename Range>
auto normalise(Range& values)
{
    return normalise(std::begin(values), std::end(values));
}

template <typename ForwardIterator>
auto normalise_logs(const ForwardIterator first, const ForwardIterator last)
{
    const auto norm = log_sum_exp(first, last);
    std::for_each(first, last, [norm] (auto& value) { value -= norm; });
    return norm;
}

template <typename Range>
auto normalise_logs(Range& logs)
{
    return normalise_logs(std::begin(logs), std::end(logs));
}

template <typename ForwardIterator>
auto normalise_exp(const ForwardIterator first, const ForwardIterator last)
{
    const auto norm = log_sum_exp(first, last);
    std::transform(first, last, first, [norm] (auto& value) { return std::exp(value - norm); });
    return norm;
}

template <typename Range>
auto normalise_exp(Range& logs)
{
    return normalise_exp(std::begin(logs), std::end(logs));
}

} // namespace maths
} // namespace octopus

#endif
