//
//  maths.hpp
//  Octopus
//
//  Created by Daniel Cooke on 23/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__maths__
#define __Octopus__maths__

#include <vector>
#include <cstddef>     // size_t
#include <cmath>       // std::exp, std::log, std::sqrt
#include <numeric>     // std::accumulate, std::iota, std::inner_product
#include <algorithm>   // std::max, std::max_element, std::transform, std::all_of
#include <type_traits> // std::enable_if, std::is_integral
#include <iterator>    // std::begin, std::end, std::cbegin, std::cend, std::distance
#include <functional>  // std::plus, std::minus

#include <iostream> // TEST

#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/special_functions/sign.hpp>
#include <boost/math/special_functions/digamma.hpp>

namespace Octopus { namespace Maths {

template <typename T, typename = typename std::enable_if_t<!std::is_integral<T>::value, bool>>
bool almost_equal(T lhs, T rhs, int ulp = 1)
{
    return lhs == rhs || std::abs(lhs - rhs) < std::numeric_limits<T>::epsilon() * std::abs(lhs + rhs) * ulp;
}

template <typename T, typename = typename std::enable_if_t<!std::is_integral<T>::value, bool>>
bool almost_zero(T x, int ulp = 1)
{
    return almost_equal(x, T {0}, ulp);
}

template <typename T, typename = typename std::enable_if_t<!std::is_integral<T>::value, bool>>
bool almost_one(T x, int ulp = 1)
{
    return almost_equal(x, T {1}, ulp);
}

template <typename InputIterator>
double mean(InputIterator first, InputIterator last)
{
    return std::accumulate(first, last, 0.0) / std::distance(first, last);
}

template <typename InputIterator>
double stdev(InputIterator first, InputIterator last)
{
    auto m = mean(first, last);
    auto n = std::distance(first, last);
    std::vector<double> diff(n);
    std::transform(first, last, std::begin(diff), std::bind2nd(std::minus<double>(), m));
    return std::sqrt(std::inner_product(diff.begin(), diff.end(), diff.begin(), 0.0) / n);
}

template <typename RealType>
inline constexpr RealType exp_maclaurin(RealType x) {
    return (6 + x * (6 + x * (3 + x))) * 0.16666666;
}

template <typename RealType>
inline constexpr RealType mercator(RealType x) {
    return x - x * x / 2 + x * x * x / 3;
}

template <typename RealType>
inline RealType log_sum_exp(RealType log_a, RealType log_b)
{
    auto r = std::minmax(log_a, log_b);
    return r.second + std::log(1 + std::exp(r.first - r.second));
}

template <typename RealType>
inline RealType log_sum_exp(RealType log_a, RealType log_b, RealType log_c)
{
    auto max = std::max({log_a, log_b, log_c});
    return max + std::log(std::exp(log_a - max) + std::exp(log_b - max) + std::exp(log_c - max));
}

template <typename RealType>
inline RealType log_sum_exp(RealType log_a, RealType log_b, RealType log_c, RealType log_d)
{
    auto max = std::max({log_a, log_b, log_c, log_d});
    return max + std::log(std::exp(log_a - max) + std::exp(log_b - max) + std::exp(log_c - max) +
                          std::exp(log_d - max));
}

template <typename RealType>
inline RealType log_sum_exp(std::initializer_list<RealType> il)
{
    auto max = std::max(il);
    RealType exp_sum {};
    for (auto x : il) {
        exp_sum += std::exp(x - max);
    }
    return max + std::log(exp_sum);
}

template <typename RealType, typename Iterator>
inline RealType log_sum_exp(Iterator first, Iterator last)
{
    auto max = *std::max_element(first, last);
    RealType exp_sum {};
    for (; first != last; ++first) {
        exp_sum += std::exp(*first - max);
    }
    return max + std::log(exp_sum);
}

template <typename RealType, typename IntegerType,
          typename = typename std::enable_if<std::is_floating_point<RealType>::value>::type,
          typename = typename std::enable_if<std::is_integral<IntegerType>::value>::type>
inline RealType log_factorial(IntegerType x)
{
    if (x == 0 || x == 1) return 0;
    if (x == 2) return std::log(2);
    if (x == 3) return std::log(6);
    if (x == 4) return std::log(24);
    
    if (x > 100) {
        return x * std::log(x) - x; // Stirling's approximation
    } else {
        std::vector<IntegerType> lx(x);
        std::iota(std::begin(lx), std::end(lx), 1);
        std::vector<RealType> tx(x);
        std::transform(std::cbegin(lx), std::cend(lx), std::begin(tx),
                       [] (IntegerType a) { return std::log(static_cast<RealType>(a)); });
        return std::accumulate(std::cbegin(tx), std::cend(tx), RealType {});
    }
}

template <typename RealType, typename InputIterator>
inline RealType log_beta(InputIterator first, InputIterator last)
{
    return std::accumulate(first, last, RealType {},
                           [] (auto v, auto x) { return v + boost::math::lgamma(x); })
            - boost::math::lgamma(std::accumulate(first, last, RealType {}));
}

template <typename RealType, typename InputIterator1, typename InputIterator2>
inline RealType log_dirichlet(InputIterator1 firstalpha, InputIterator1 lastalpha, InputIterator2 firstpi)
{
    return std::inner_product(firstalpha, lastalpha, firstpi, RealType {}, std::plus<RealType>(),
                              [] (auto a, auto p) { return (a - 1) * std::log(p); })
            - log_beta<RealType>(firstalpha, lastalpha);
}

template <typename RealType, typename IntegerType>
inline RealType log_multinomial_coefficient(std::initializer_list<IntegerType> il)
{
    using std::begin; using std::end; using std::cbegin; using std::cend; using std::accumulate;
    std::vector<RealType> denoms(il.size());
    std::transform(cbegin(il), cend(il), begin(denoms), log_factorial<RealType, IntegerType>);
    return log_factorial<RealType>(accumulate(cbegin(il), cend(il), 0))
            - accumulate(cbegin(denoms), cend(denoms), RealType {});
}

template <typename RealType, typename Iterator>
inline RealType log_multinomial_coefficient(Iterator first, Iterator last)
{
    using IntegerType = typename Iterator::value_type;
    std::vector<RealType> denoms(std::distance(first, last));
    std::transform(first, last, std::begin(denoms), log_factorial<RealType, IntegerType>);
    return log_factorial<RealType, IntegerType>(std::accumulate(first, last, IntegerType {})) -
                std::accumulate(denoms.cbegin(), denoms.cend(), RealType {});
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

template <typename IntegerType, typename RealType>
inline RealType multinomial_pdf(const std::vector<IntegerType>& z, const std::vector<RealType>& p)
{
    RealType r {1};
    
    for (size_t i {}; i < z.size(); ++i) {
        r *= std::pow(p[i], z[i]);
    }
    
    return multinomial_coefficient<IntegerType, RealType>(std::cbegin(z), std::cend(z)) * r;
}

// Returns approximate y such that digamma(y) = x
template <typename RealType>
inline RealType digamma_inv(RealType x)
{
    RealType l {1.0};
    auto y = std::exp(x);
    
    while (l > 10e-8) {
        y += l * boost::math::sign(x - boost::math::digamma<RealType>(y));
        l /= 2;
    }
    
    return y;
}

template <typename RealType>
RealType dirichlet_multinomial(RealType z1, RealType z2, RealType a1, RealType a2)
{
    auto z_0 = z1 + z2;
    auto a_0 = a1 + a2;
    auto z_m = boost::math::factorial<RealType>(z1) * boost::math::factorial<RealType>(z2);
    
    return (boost::math::factorial<RealType>(z_0) / z_m) *
            (boost::math::tgamma(a_0) / boost::math::tgamma(z_0 + a_0)) *
            (boost::math::tgamma<RealType>(z1 + a1) * boost::math::tgamma<RealType>(z2 + a2)) /
            (boost::math::tgamma<RealType>(a1) + boost::math::tgamma<RealType>(a2));
}

template <typename RealType>
RealType dirichlet_multinomial(RealType z1, RealType z2, RealType z3, RealType a1, RealType a2, RealType a3)
{
    auto z_0 = z1 + z2 + z3;
    auto a_0 = a1 + a2 + a3;
    auto z_m = boost::math::factorial<RealType>(z1) * boost::math::factorial<RealType>(z2) *
                boost::math::factorial<RealType>(z3);
    
    return (boost::math::factorial<RealType>(z_0) / z_m) *
            (boost::math::tgamma(a_0) / boost::math::tgamma(z_0 + a_0)) *
            (boost::math::tgamma<RealType>(z1 + a1) * boost::math::tgamma<RealType>(z2 + a2) *
             boost::math::tgamma<RealType>(z3 + a3)) /
            (boost::math::tgamma<RealType>(a1) + boost::math::tgamma<RealType>(a2) + boost::math::tgamma<RealType>(a3));
}

template <typename RealType>
RealType dirichlet_multinomial(const std::vector<RealType>& z, const std::vector<RealType>& a)
{
    auto z_0 = std::accumulate(std::cbegin(z), std::cend(z), RealType {});
    auto a_0 = std::accumulate(std::cbegin(a), std::cend(a), RealType {});
    
    RealType z_m {1};
    for (auto z_i : z) {
        z_m *= boost::math::factorial<RealType>(z_i);
    }
    
    RealType g {1};
    for (size_t i {0}; i < z.size(); ++i) {
        g *= boost::math::tgamma<RealType>(z[i] + a[i]) / boost::math::tgamma<RealType>(a[i]);
    }
    
    return (boost::math::factorial<RealType>(z_0) / z_m) *
            (boost::math::tgamma(a_0) / boost::math::tgamma(z_0 + a_0)) * g;
}

template <typename RealType>
RealType beta_binomial(RealType k, RealType n, RealType alpha, RealType beta)
{
    return dirichlet_multinomial<RealType>(k, n - k, alpha, beta);
}

namespace detail {
template <typename RealType>
bool is_mldp_converged(std::vector<RealType>& lhs, const std::vector<RealType>& rhs, const RealType epsilon = 0.0001)
{
    using std::cbegin; using std::cend; using std::begin;
    std::transform(cbegin(lhs), cend(lhs), cbegin(rhs), begin(lhs), [] (auto a, auto b) { return std::abs(a - b); });
    return std::all_of(cbegin(lhs), cend(lhs), [epsilon] (auto x) { return x < epsilon; });
}
} // namespace detail


template <typename RealType>
std::vector<RealType>
maximum_likelihood_dirichlet_params(std::vector<RealType> pi, const RealType precision,
                                    const unsigned max_iterations = 100, const RealType epsilon = 0.0001)
{
    using std::cbegin; using std::cend; using std::begin;
    
    std::transform(cbegin(pi), cend(pi), begin(pi), [] (auto p) { return std::log(p); });
    
    const auto l = pi.size();
    std::vector<RealType> result(l, 1.0 / l), curr_result(l, 1.0 / l), means(l, 1.0 / l);
    
    for (unsigned n {}; n < max_iterations; ++n) {
        RealType v {};
        
        for (size_t j {}; j < l; ++j) {
            v += means[j] * (pi[j] - boost::math::digamma<RealType>(precision * means[j]));
        }
        
        for (size_t k {}; k < l; ++k) {
            curr_result[k] = digamma_inv<RealType>(pi[k] - v);
            means[k]       = curr_result[k] / std::accumulate(cbegin(curr_result), cend(curr_result), RealType {});
        }
        
        if (detail::is_mldp_converged(result, curr_result, epsilon)) return curr_result;
        
        result = curr_result;
    }
    
    return result;
}

template <typename MapType>
inline
typename MapType::key_type
sum_keys(const MapType& map)
{
    return std::accumulate(std::cbegin(map), std::cend(map), typename MapType::key_type {},
                           [] (const auto previous, const auto& p) { return previous + p.first; });
}

template <typename ResultType, typename MapType, typename UnaryOperation>
inline
ResultType
sum_keys(const MapType& map, UnaryOperation op)
{
    return std::accumulate(std::cbegin(map), std::cend(map), ResultType {},
                           [op] (const auto previous, const auto& p) { return previous + op(p.first); });
}

template <typename MapType>
inline
typename MapType::mapped_type
sum_values(const MapType& map)
{
    return std::accumulate(std::cbegin(map), std::cend(map), typename MapType::mapped_type {},
                           [] (const auto previous, const auto& p) { return previous + p.second; });
}

template <typename ResultType, typename MapType, typename UnaryOperation>
inline
ResultType
sum_values(const MapType& map, UnaryOperation op)
{
    return std::accumulate(std::cbegin(map), std::cend(map), ResultType {},
                           [op] (const auto previous, const auto& p) { return previous + op(p.second); });
}

template <typename Map>
inline
size_t sum_sizes(const Map& map)
{
    return std::accumulate(std::cbegin(map), std::cend(map), size_t {},
                           [] (const auto& p, const auto& v) { return p + v.second.size(); });
}

} // namespace Maths
} // namespace Octopus

#endif /* defined(__Octopus__maths__) */
