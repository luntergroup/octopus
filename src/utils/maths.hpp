// Copyright (c) 2016 Daniel Cooke
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

namespace octopus { namespace maths
{
namespace constants
{
    template <typename T = double>
    constexpr T ln10Div10 = T {0.230258509299404568401799145468436420760110148862877297603};
}

template <typename RealType,
          typename = std::enable_if_t<std::is_floating_point<RealType>::value>>
RealType round(const RealType val, const unsigned precision = 2)
{
    const auto factor = std::pow(RealType {10.0}, precision);
    return std::round(val * factor) / factor;
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


template <typename RealType>
constexpr RealType exp_maclaurin(const RealType x) {
    return (6 + x * (6 + x * (3 + x))) * 0.16666666;
}

template <typename RealType>
constexpr RealType mercator(const RealType x) {
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

template <typename InputIt, typename UnaryOperation>
double stdev(InputIt first, InputIt last, UnaryOperation unary_op)
{
    const auto m = mean(first, last, unary_op);
    const auto n = std::distance(first, last);
    std::vector<double> diff(n);
    std::transform(first, last, std::begin(diff), [&] (const auto& x) { return unary_op(x) - m; });
    return std::sqrt(std::inner_product(std::begin(diff), std::end(diff), std::begin(diff), 0.0) / n);
}

template <typename InputIt>
auto stdev(InputIt first, InputIt last)
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

template <typename RealType, typename InputIt>
RealType rmq(InputIt first, InputIt last)
{
    return std::sqrt((std::inner_product(first, last, first, RealType {0}))
                     / static_cast<RealType>(std::distance(first, last)));
}

template <typename RealType, typename Container>
RealType rmq(const Container& values)
{
    return rmq<RealType>(std::cbegin(values), std::cend(values));
}

template <typename RealType,
          typename = std::enable_if_t<std::is_floating_point<RealType>::value>>
RealType log_sum_exp(const RealType a, const RealType b)
{
    const auto r = std::minmax(a, b);
    return r.second + std::log(1.0 + std::exp(r.first - r.second));
}

template <typename RealType,
          typename = std::enable_if_t<std::is_floating_point<RealType>::value>>
RealType log_sum_exp(const RealType a, const RealType b, const RealType c)
{
    const auto max = std::max({a, b, c});
    return max + std::log(std::exp(a - max) + std::exp(b - max) + std::exp(c - max));
}

template <typename RealType,
          typename = std::enable_if_t<std::is_floating_point<RealType>::value>>
RealType log_sum_exp(std::initializer_list<RealType> il)
{
    const auto max = std::max(il);
    return max + std::log(std::accumulate(std::cbegin(il), std::cend(il), RealType {0},
                                          [max] (const auto curr, const auto x) {
                                              return curr + std::exp(x - max);
                                          }));
}

template <typename ForwardIt,
          typename = std::enable_if_t<!std::is_floating_point<ForwardIt>::value>>
auto log_sum_exp(ForwardIt first, ForwardIt last)
{
    assert(first != last);
    using RealType = typename std::iterator_traits<ForwardIt>::value_type;
    const auto max = *std::max_element(first, last);
    return max + std::log(std::accumulate(first, last, RealType {0},
                                          [max] (const auto curr, const auto x) {
                                              return curr + std::exp(x - max);
                                          }));
}

template <typename Container>
auto log_sum_exp(const Container& values)
{
    return log_sum_exp(std::cbegin(values), std::cend(values));
}

template <typename RealType, typename IntegerType,
          typename = std::enable_if_t<std::is_floating_point<RealType>::value>,
          typename = std::enable_if_t<std::is_integral<IntegerType>::value>>
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
        std::iota(std::begin(lx), std::end(lx), 1);
        std::vector<RealType> tx(x);
        std::transform(std::cbegin(lx), std::cend(lx), std::begin(tx),
                       [] (IntegerType a) { return std::log(static_cast<RealType>(a)); });
        return std::accumulate(std::cbegin(tx), std::cend(tx), RealType {0});
    }
}

template <typename ForwardIt>
auto log_beta(const ForwardIt first, const ForwardIt last)
{
    using T = std::decay_t<typename std::iterator_traits<ForwardIt>::value_type>;
    static_assert(std::is_floating_point<T>::value,
                  "log_beta is only defined for floating point types.");
    return std::accumulate(first, last, T {0},
                           [] (const auto curr, const auto x) {
                               return curr + boost::math::lgamma(x);
                           }) - boost::math::lgamma(std::accumulate(first, last, T {0}));
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
    return std::inner_product(firstalpha, lastalpha, firstpi, T {0}, std::plus<void> {},
                              [] (const auto a, const auto p) { return (a - 1) * std::log(p); })
            - log_beta(firstalpha, lastalpha);
}

template <typename Container1, typename Container2>
auto log_dirichlet(const Container1& alpha, const Container2& pi)
{
    return log_dirichlet(std::cbegin(alpha), std::cend(alpha), std::cbegin(pi));
}

template <typename RealType, typename IntegerType>
inline RealType log_multinomial_coefficient(std::initializer_list<IntegerType> il)
{
    using std::begin; using std::end; using std::cbegin; using std::cend; using std::accumulate;
    std::vector<RealType> denoms(il.size());
    std::transform(cbegin(il), cend(il), begin(denoms), log_factorial<RealType, IntegerType>);
    return log_factorial<RealType>(accumulate(cbegin(il), cend(il), 0))
            - accumulate(cbegin(denoms), cend(denoms), RealType {0});
}

template <typename RealType, typename Iterator>
inline RealType log_multinomial_coefficient(Iterator first, Iterator last)
{
    using IntegerType = typename Iterator::value_type;
    std::vector<RealType> denoms(std::distance(first, last));
    std::transform(first, last, std::begin(denoms), log_factorial<RealType, IntegerType>);
    return log_factorial<RealType, IntegerType>(std::accumulate(first, last, IntegerType {0}))
            - std::accumulate(std::cbegin(denoms), std::cend(denoms), RealType {0});
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

template <typename IntegerType, typename RealType>
inline RealType multinomial_pdf(const std::vector<IntegerType>& z, const std::vector<RealType>& p)
{
    RealType r {1};
    for (std::size_t i {0}; i < z.size(); ++i) {
        r *= std::pow(p[i], z[i]);
    }
    return multinomial_coefficient<IntegerType, RealType>(std::cbegin(z), std::cend(z)) * r;
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

template <typename RealType>
RealType dirichlet_multinomial(const RealType z1, const RealType z2, const RealType a1, const RealType a2)
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
RealType dirichlet_multinomial(const RealType z1, const RealType z2, const RealType z3,
                               const RealType a1, const RealType a2, const RealType a3)
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
    auto z_0 = std::accumulate(std::cbegin(z), std::cend(z), RealType {0});
    auto a_0 = std::accumulate(std::cbegin(a), std::cend(a), RealType {0});
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
}

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
            means[k]       = curr_result[k] / std::accumulate(std::cbegin(curr_result),
                                                              std::cend(curr_result),
                                                              RealType {0});
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
NumericType probability_to_phred(const RealType p)
{
    return static_cast<NumericType>(-10.0 * std::log10(std::max(1.0 - p, std::numeric_limits<RealType>::epsilon())));
}

template <typename NumericType = float, typename RealType,
          typename = std::enable_if_t<std::is_floating_point<RealType>::value>>
NumericType probability_to_phred(const RealType p, const unsigned precision)
{
    return round(static_cast<NumericType>(-10.0 * std::log10(std::max(1.0 - p, std::numeric_limits<RealType>::epsilon()))), precision);
}

template <typename RealType = double, typename NumericType>
RealType phred_to_probability(const NumericType phred)
{
    return 1.0 - std::pow(10.0, -1.0 * static_cast<RealType>(phred) / 10.0);
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
RealType beta_cdf_complement(const RealType a, const RealType b, const RealType x)
{
    const boost::math::beta_distribution<> beta_dist {a, b};
    return boost::math::cdf(boost::math::complement(beta_dist, x));
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
    return std::make_pair(boost::math::ibeta_inv(a, b, c),
                          boost::math::ibeta_inv(a, b, c + mass));
}

} // namespace detail

template <typename RealType>
std::pair<RealType, RealType>
beta_hdi(RealType a, RealType b, const RealType mass = 0.99)
{
    static_assert(std::is_floating_point<RealType>::value,
                  "beta_hdi only works for floating point types");
    
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

template <typename Container>
void log_each(Container& values)
{
    for (auto& v : values) v = std::log(v);
}

template <typename Container>
void exp_each(Container& values)
{
    for (auto& v : values) v = std::exp(v);
}

template <typename Container>
auto normalise_logs(Container& logs)
{
    const auto norm = log_sum_exp(logs);
    for (auto& p : logs) p -= norm;
    return norm;
}

template <typename Container>
auto normalise_exp(Container& logs)
{
    const auto norm = log_sum_exp(logs);
    for (auto& p : logs) p = std::exp(p -= norm);
    return norm;
}
} // namespace maths
} // namespace octopus

#endif
