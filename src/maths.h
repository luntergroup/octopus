//
//  maths.h
//  Octopus
//
//  Created by Daniel Cooke on 23/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__maths__
#define __Octopus__maths__

#include <vector>
#include <cstddef>     // std::size_t
#include <cmath>       // std::exp, std::log
#include <numeric>     // std::accumulate, std::iota
#include <algorithm>   // std::max, std::max_element, std::transform
#include <type_traits> // std::enable_if, std::is_integral

#include <boost/math/special_functions/gamma.hpp>
#include <boost/math/special_functions/factorials.hpp>
#include <boost/math/special_functions/sign.hpp>
#include <boost/math/special_functions/digamma.hpp>

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
    for (const auto& x : il) {
        exp_sum += std::exp(x - max);
    }
    return max + std::log(exp_sum);
}

template <typename RealType, typename Iterator>
inline RealType log_sum_exp(Iterator first, Iterator last)
{
    auto max = *std::max_element(first, last);
    RealType exp_sum {};
    while (first != last) {
        exp_sum += std::exp(*first - max);
        ++first;
    }
    return max + std::log(exp_sum);
}

template <typename RealType, typename IntegerType,
          typename = typename std::enable_if<std::is_floating_point<RealType>::value>::type,
          typename = typename std::enable_if<std::is_integral<IntegerType>::value>::type>
inline
RealType log_factorial(IntegerType x)
{
    if (x == 0 || x == 1) return 0;
    if (x == 2) return std::log(2);
    if (x == 3) return std::log(6);
    if (x == 4) return std::log(24);
    
    if (x > 100) {
        return x * std::log(x) - x; // Stirling's approximation
    } else {
        std::vector<IntegerType> lx(x);
        std::iota(lx.begin(), lx.end(), 1);
        std::vector<RealType> tx(x);
        std::transform(lx.cbegin(), lx.cend(), tx.begin(),
                       [] (IntegerType a) { return std::log(static_cast<RealType>(a)); });
        return std::accumulate(tx.cbegin(), tx.cend(), static_cast<RealType>(0));
    }
}

template <typename RealType, typename IntegerType>
inline
RealType
log_multinomial_coefficient(std::initializer_list<IntegerType> il)
{
    std::vector<RealType> denoms(il.size());
    std::transform(il.begin(), il.end(), denoms.begin(), log_factorial<RealType, IntegerType>);
    return log_factorial<RealType>(std::accumulate(il.begin(), il.end(), 0)) -
            std::accumulate(denoms.cbegin(), denoms.cend(), static_cast<RealType>(0));
}

template <typename RealType, typename Iterator>
inline
RealType
log_multinomial_coefficient(Iterator first, Iterator last)
{
    std::vector<RealType> denoms(std::distance(first, last));
    using IntegerType = typename Iterator::value_type;
    std::transform(first, last, denoms.begin(), log_factorial<RealType, IntegerType>);
    return log_factorial<RealType>(std::accumulate(first, last, 0)) -
            std::accumulate(denoms.cbegin(), denoms.cend(), static_cast<RealType>(0));
}

template <typename RealType, typename IntegerType>
inline
IntegerType
multinomial_coefficient(std::initializer_list<IntegerType> il)
{
    return static_cast<IntegerType>(std::exp(log_multinomial_coefficient<RealType, IntegerType>(std::move(il))));
}

template <typename IntegerType, typename RealType, typename Iterator>
inline
IntegerType
multinomial_coefficient(Iterator first, Iterator last)
{
    return static_cast<IntegerType>(std::exp(log_multinomial_coefficient<RealType, IntegerType>(first, last)));
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
    auto z_0 = std::accumulate(z.cbegin(), z.cend(), RealType {});
    auto a_0 = std::accumulate(a.cbegin(), a.cend(), RealType {});
    
    RealType z_m {1};
    for (auto z_i : z) {
        z_m *= boost::math::factorial<RealType>(z_i);
    }
    
    RealType g {1};
    for (std::size_t i {0}; i < z.size(); ++i) {
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

// Returns approximate y such that digamma(y) = x
template <typename RealType>
RealType digamma_inv(RealType x)
{
    RealType l {1.0};
    auto y = std::exp(x);
    
    while (l > 10e-8) {
        y += l * boost::math::sign(x - boost::math::digamma<RealType>(y));
        l /= 2;
    }
    
    return y;
}

template <typename MapType>
inline
typename MapType::mapped_type
sum_values(const MapType& map)
{
    typename MapType::mapped_type result {};
    
    for (const auto& map_pair : map) {
        result += map_pair.second;
    }
    
    return result;
}

template <typename ResultType, typename MapType, typename UnaryOperation>
inline
ResultType
sum_values(const MapType& map, UnaryOperation op)
{
    ResultType result {};
    
    for (const auto& map_pair : map) {
        result += op(map_pair.second);
    }
    
    return result;
}

#endif /* defined(__Octopus__maths__) */
