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
#include <cmath>
#include <numeric>   // std::accumulate
#include <algorithm> // std::max, std::max_element

template <typename T>
inline constexpr T factorial(T x)
{
    return (x == 0) ? 1 : x * factorial(x - 1);
}

template <typename T>
inline T multinomial_coefficient(const std::vector<T>& xs)
{
    auto x_0 = std::accumulate(xs.cbegin(), xs.cend(), T {});
    return x_0;
}

template <typename T>
inline constexpr T exp_maclaurin(T x) {
    return (6 + x * (6 + x * (3 + x))) * 0.16666666;
}

template <typename T>
inline constexpr T mercator(T x) {
    return x - x * x / 2 + x * x * x / 3;
}

template <typename T>
inline T log_sum_exp(T log_a, T log_b)
{
    auto r = std::minmax(log_a, log_b);
    return r.second + std::log(1 + std::exp(r.first - r.second));
}

template <typename T>
inline T log_sum_exp(T log_a, T log_b, T log_c)
{
    auto max = std::max({log_a, log_b, log_c});
    return max + std::log(std::exp(log_a - max) + std::exp(log_b - max) + std::exp(log_c - max));
}

template <typename T>
inline T log_sum_exp(T log_a, T log_b, T log_c, T log_d)
{
    auto max = std::max({log_a, log_b, log_c, log_d});
    return max + std::log(std::exp(log_a - max) + std::exp(log_b - max) + std::exp(log_c - max) +
                          std::exp(log_d - max));
}

template <typename T>
inline T log_sum_exp(std::initializer_list<T> il)
{
    auto max = std::max(il);
    T exp_sum {};
    for (const auto& x : il) {
        exp_sum += std::exp(x - max);
    }
    return max + std::log(exp_sum);
}

template <typename T, typename Iterator>
inline T log_sum_exp(Iterator first, Iterator last)
{
    auto max = *std::max_element(first, last);
    T exp_sum {};
    while (first != last) {
        exp_sum += std::exp(*first - max);
        ++first;
    }
    return max + std::log(exp_sum);
}

//template <typename RealType>
//RealType dirichlet(const std::vector<RealType>& x, const std::vector<RealType>& a)
//{
//    
//}
//
//template <typename RealType>
//RealType dirichlet_multinomial(const std::vector<RealType>& x, const std::vector<RealType>& a)
//{
//    auto a_0 = std::accumulate(a.cbegin(), a.cend(), RealType{});
//    auto s_0 = std::accumulate(x.cbegin(), x.cend(), RealType{});
//    
//    RealType c {};
//    for (const auto x_k : x) {
//        c *= factorial(x_k);
//    }
//    
//    RealType 
//}

#endif /* defined(__Octopus__maths__) */
