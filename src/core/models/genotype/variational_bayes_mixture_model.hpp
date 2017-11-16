// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef variational_bayes_mixture_model_hpp
#define variational_bayes_mixture_model_hpp

#include <array>
#include <vector>
#include <algorithm>
#include <numeric>
#include <iterator>
#include <cstddef>
#include <utility>
#include <cassert>

#include <boost/math/special_functions/digamma.hpp>

#include "utils/maths.hpp"
#include "core/models/haplotype_likelihood_cache.hpp"

/**
 *
 * This file contains an implementation of the Variational Bayes mixture model
 * used by some genotype models. The notation follows the documentation.
 *
 */

namespace octopus { namespace model {

// Types needed for Variational Bayes model

struct VariationalBayesParameters
{
    double epsilon;
    unsigned max_iterations;
};

using ProbabilityVector    = std::vector<double>;
using LogProbabilityVector = std::vector<double>;

template <std::size_t K>
using VBAlpha = std::array<double, K>;
template <std::size_t K>
using VBAlphaVector = std::vector<VBAlpha<K>>;

class VBReadLikelihoodArray
{
public:
    using BaseType = HaplotypeLikelihoodCache::LikelihoodVector;
    
    VBReadLikelihoodArray() = default;
    
    explicit VBReadLikelihoodArray(const BaseType&);
    
    VBReadLikelihoodArray(const VBReadLikelihoodArray&)            = default;
    VBReadLikelihoodArray& operator=(const VBReadLikelihoodArray&) = default;
    VBReadLikelihoodArray(VBReadLikelihoodArray&&)                 = default;
    VBReadLikelihoodArray& operator=(VBReadLikelihoodArray&&)      = default;
    
    ~VBReadLikelihoodArray() = default;
    
    void operator=(const BaseType&);
    void operator=(std::reference_wrapper<const BaseType>);
    std::size_t size() const noexcept;
    BaseType::const_iterator begin() const noexcept;
    BaseType::const_iterator end() const noexcept;
    double operator[](const std::size_t n) const noexcept;

private:
    const BaseType* likelihoods;
};

template <std::size_t K>
using VBGenotype = std::array<VBReadLikelihoodArray, K>;
template <std::size_t K>
using VBGenotypeVector = std::vector<VBGenotype<K>>;
template <std::size_t K>
using VBReadLikelihoodMatrix = std::vector<VBGenotypeVector<K>>;

template <std::size_t K>
using VBTau = std::array<double, K>;
template <std::size_t K>
using VBResponsabilityVector = std::vector<VBTau<K>>;
template <std::size_t K>
using VBResponsabilityMatrix = std::vector<VBResponsabilityVector<K>>;

template <std::size_t K>
struct VBLatents
{
    ProbabilityVector genotype_posteriors;
    LogProbabilityVector genotype_log_posteriors;
    VBAlphaVector<K> alphas;
    VBResponsabilityMatrix<K> responsabilities;
};

// Main VB mathod

template <std::size_t K>
std::pair<VBLatents<K>, double>
run_variational_bayes(const VBAlphaVector<K>& prior_alphas,
                      const LogProbabilityVector& genotype_log_priors,
                      const VBReadLikelihoodMatrix<K>& log_likelihoods,
                      const VariationalBayesParameters& params,
                      std::vector<LogProbabilityVector> seeds);

namespace detail {

inline ProbabilityVector& exp(const LogProbabilityVector& log_probabilities, ProbabilityVector& result) noexcept
{
    std::transform(std::cbegin(log_probabilities), std::cend(log_probabilities), std::begin(result),
                   [] (const auto lp) noexcept { return std::exp(lp); });
    return result;
}

inline ProbabilityVector exp(const LogProbabilityVector& log_probabilities)
{
    ProbabilityVector result(log_probabilities.size());
    return exp(log_probabilities, result);
}

inline auto sum(const VBAlpha<2>& alpha) noexcept
{
    return alpha[0] + alpha[1];
}

inline auto sum(const VBAlpha<3>& alpha) noexcept
{
    return alpha[0] + alpha[1] + alpha[2];
}

template <std::size_t K>
auto sum(const VBAlpha<K>& alpha) noexcept
{
    return std::accumulate(std::cbegin(alpha), std::cend(alpha), 0.0);
}

inline auto digamma_diff(const double a, const double b)
{
    using boost::math::digamma;
    return digamma(a) - digamma(b);
}

inline double log_sum_exp(const std::array<double, 1>& logs)
{
    return logs[0];
}

inline double log_sum_exp(const std::array<double, 2>& logs)
{
    return maths::log_sum_exp(logs[0], logs[1]);
}

inline double log_sum_exp(const std::array<double, 3>& logs)
{
    return maths::log_sum_exp(logs[0], logs[1], logs[2]);
}

template <std::size_t K>
double log_sum_exp(const std::array<double, K>& logs)
{
    return maths::log_sum_exp(logs);
}

template <std::size_t K>
double expectation(const ProbabilityVector& distribution, const VBGenotypeVector<K>& likelihoods,
                   const unsigned k, const std::size_t n) noexcept
{
    return std::inner_product(std::cbegin(distribution), std::cend(distribution),
                              std::cbegin(likelihoods), 0.0, std::plus<> {},
                              [k, n](const auto p, const auto& haplotype_likelihoods) noexcept {
                                  return p * haplotype_likelihoods[k][n];
                              });
}

template <std::size_t K>
VBResponsabilityVector<K>
init_responsabilities(const VBAlpha<K>& prior_alphas,
                      const ProbabilityVector& genotype_probabilities,
                      const VBGenotypeVector<K>& read_likelihoods)
{
    std::array<double, K> al; // no need to keep recomputing this
    const auto a0 = sum(prior_alphas);
    for (unsigned k {0}; k < K; ++k) {
        al[k] = digamma_diff(prior_alphas[k], a0);
    }
    const auto N = read_likelihoods[0][0].size(); // num reads in sample s
    VBResponsabilityVector<K> result(N);
    std::array<double, K> ln_rho;
    for (std::size_t n {0}; n < N; ++n) {
        for (unsigned k {0}; k < K; ++k) {
            ln_rho[k] = al[k] + expectation(genotype_probabilities, read_likelihoods, k, n);
        }
        const auto ln_rho_norm = log_sum_exp(ln_rho);
        for (unsigned k {0}; k < K; ++k) {
            result[n][k] = std::exp(ln_rho[k] - ln_rho_norm);
        }
    }
    return result;
}

template <std::size_t K>
VBResponsabilityMatrix<K>
init_responsabilities(const VBAlphaVector<K>& prior_alphas,
                      const ProbabilityVector& genotype_probabilities,
                      const VBReadLikelihoodMatrix<K>& read_likelihoods)
{
    assert(!read_likelihoods.empty());
    assert(prior_alphas.size() == read_likelihoods.size());
    const auto S = read_likelihoods.size(); // num samples
    VBResponsabilityMatrix<K> result{};
    result.reserve(S);
    for (std::size_t s {0}; s < S; ++s) {
        result.push_back(init_responsabilities(prior_alphas[s], genotype_probabilities, read_likelihoods[s]));
    }
    return result;
}

template <std::size_t K>
void update_responsabilities(VBResponsabilityVector<K>& result,
                             const VBAlpha<K>& posterior_alphas,
                             const ProbabilityVector& genotype_probabilities,
                             const VBGenotypeVector<K>& read_likelihoods)
{
    std::array<double, K> al;
    const auto a0 = sum(posterior_alphas);
    for (unsigned k {0}; k < K; ++k) {
        al[k] = digamma_diff(posterior_alphas[k], a0);
    }
    const auto N = read_likelihoods[0][0].size();
    std::array<double, K> ln_rho;
    for (std::size_t n {0}; n < N; ++n) {
        for (unsigned k {0}; k < K; ++k) {
            ln_rho[k] = al[k] + expectation(genotype_probabilities, read_likelihoods, k, n);
        }
        const auto ln_rho_norm = log_sum_exp(ln_rho);
        for (unsigned k {0}; k < K; ++k) {
            result[n][k] = std::exp(ln_rho[k] - ln_rho_norm);
        }
    }
}

// same as init_responsabilities but in-place
template <std::size_t K>
void update_responsabilities(VBResponsabilityMatrix<K>& result,
                             const VBAlphaVector<K>& posterior_alphas,
                             const ProbabilityVector& genotype_probabilities,
                             const VBReadLikelihoodMatrix<K>& read_likelihoods)
{
    const auto S = read_likelihoods.size();
    for (std::size_t s {0}; s < S; ++s) {
        update_responsabilities(result[s], posterior_alphas[s], genotype_probabilities, read_likelihoods[s]);
    }
}

template <std::size_t K>
auto sum(const VBResponsabilityVector<K>& taus, const unsigned k) noexcept
{
    return std::accumulate(std::cbegin(taus), std::cend(taus), 0.0,
                           [k](const auto curr, const auto& tau) noexcept {
                               return curr + tau[k];
                           });
}

template <std::size_t K>
void update_alpha(VBAlpha<K>& alpha, const VBAlpha<K>& prior_alpha,
                  const VBResponsabilityVector<K>& taus) noexcept
{
    for (unsigned k {0}; k < K; ++k) {
        alpha[k] = prior_alpha[k] + sum(taus, k);
    }
}

template <std::size_t K>
void update_alphas(VBAlphaVector<K>& alphas, const VBAlphaVector<K>& prior_alphas,
                   const VBResponsabilityMatrix<K>& responsabilities) noexcept
{
    const auto S = alphas.size();
    assert(S == prior_alphas.size() && S == responsabilities.size());
    for (std::size_t s {0}; s < S; ++s) {
        update_alpha(alphas[s], prior_alphas[s], responsabilities[s]);
    }
}

template <std::size_t K>
auto marginalise(const VBResponsabilityMatrix<K>& responsabilities,
                 const VBReadLikelihoodMatrix<K>& read_likelihoods,
                 const std::size_t g)
{
    double result {0};
    const auto S = read_likelihoods.size(); // num samples
    assert(S == responsabilities.size());
    for (std::size_t s {0}; s < S; ++s) {
        const auto N = read_likelihoods[s][0][0].size(); // num reads in sample s
        assert(responsabilities[s].size() == N);
        assert(responsabilities[s][0].size() == K);
        assert(read_likelihoods[s][g].size() == K);
        for (unsigned k {0}; k < K; ++k) {
            for (std::size_t n {0}; n < N; ++n) {
                result += responsabilities[s][n][k] * read_likelihoods[s][g][k][n];
            }
        }
    }
    return result;
}

template <std::size_t K>
void update_genotype_log_posteriors(LogProbabilityVector& result,
                                    const LogProbabilityVector& genotype_log_priors,
                                    const VBResponsabilityMatrix<K>& responsabilities,
                                    const VBReadLikelihoodMatrix<K>& read_likelihoods)
{
    const auto G = result.size();
    for (std::size_t g {0}; g < G; ++g) {
        result[g] = genotype_log_priors[g] + marginalise(responsabilities, read_likelihoods, g);
    }
    maths::normalise_logs(result);
}

inline auto max_change(const VBAlpha<2>& lhs, const VBAlpha<2>& rhs) noexcept
{
    return std::max(std::abs(lhs.front() - rhs.front()), std::abs(lhs.back() - rhs.back()));
}

inline auto max_change(const VBAlpha<3>& lhs, const VBAlpha<3>& rhs) noexcept
{
    return std::max({std::abs(lhs[0] - rhs[0]), std::abs(lhs[1] - rhs[1]), std::abs(lhs[2] - rhs[2])});
}

template <std::size_t K>
auto max_change(const VBAlpha<K>& lhs, const VBAlpha<K>& rhs) noexcept
{
    double result {0};
    for (std::size_t k {0}; k < K; ++k) {
        const auto curr = std::abs(lhs[k] - rhs[k]);
        if (curr > result) result = curr;
    }
    return result;
}

template <std::size_t K>
auto max_change(const VBAlphaVector<K>& prior_alphas, const VBAlphaVector<K>& posterior_alphas) noexcept
{
    const auto S = prior_alphas.size();
    assert(S == posterior_alphas.size());
    double result {0};
    for (std::size_t s {0}; s < S; ++s) {
        const auto curr = max_change(prior_alphas[s], posterior_alphas[s]);
        if (curr > result) result = curr;
    }
    return result;
}

template <std::size_t K>
std::pair<bool, double> check_convergence(const VBAlphaVector<K>& prior_alphas,
                                          const VBAlphaVector<K>& posterior_alphas,
                                          const double prev_max_change,
                                          const double epsilon) noexcept
{
    const auto new_max_change = max_change(prior_alphas, posterior_alphas);
    return std::make_pair(std::abs(new_max_change - prev_max_change) < epsilon, new_max_change);
}

// lower-bound calculation

inline auto expectation(const ProbabilityVector& genotype_posteriors,
                        const LogProbabilityVector& genotype_log_priors) noexcept
{
    return std::inner_product(std::cbegin(genotype_posteriors), std::cend(genotype_posteriors),
                              std::cbegin(genotype_log_priors), 0.0);
}

template <std::size_t K>
auto dirichlet_expectation(const VBAlpha<K>& priors, const VBAlpha<K>& posteriors)
{
    using boost::math::digamma;
    const auto da0 = digamma(sum(posteriors));
    return std::inner_product(std::cbegin(priors), std::cend(priors),
                              std::cbegin(posteriors), 0.0, std::plus<> {},
                              [da0](const auto& prior, const auto& post) {
                                  return (prior - 1) * (digamma(post) - da0);
                              }) - maths::log_beta(priors);
}

template <std::size_t K>
auto expectation(const VBAlphaVector<K>& priors, const VBAlphaVector<K>& posteriors)
{
    return std::inner_product(std::cbegin(priors), std::cend(priors),
                              std::cbegin(posteriors), 0.0, std::plus<> {},
                              [](const auto& prior, const auto& post) {
                                  return dirichlet_expectation(prior, post);
                              });
}

// E[ln p(Z_s | pi_s)]
template <std::size_t K>
auto expectation(const VBResponsabilityVector<K>& taus, const VBAlpha<K>& alpha)
{
    using boost::math::digamma;
    const auto das = digamma(sum(alpha));
    double result {0};
    for (unsigned k {0}; k < K; ++k) {
        result += (digamma(alpha[k]) - das) * sum(taus, k);
    }
    return result;
}

// sum s E[ln p(Z_s | pi_s)]
template <std::size_t K>
auto expectation(const VBResponsabilityMatrix<K>& taus, const VBAlphaVector<K>& alphas)
{
    return std::inner_product(std::cbegin(taus), std::cend(taus), std::cbegin(alphas), 0.0, std::plus<> {},
                              [](const auto& tau, const auto& alpha) {
                                  return expectation(tau, alpha);
                              });
}

template <std::size_t K>
auto expectation(const VBResponsabilityMatrix<K>& taus,
                   const VBReadLikelihoodMatrix<K>& log_likelihoods,
                   const std::size_t g)
{
    double result {0};
    for (std::size_t s {0}; s < taus.size(); ++s) {
        for (std::size_t n {0}; n < taus[s].size(); ++n) {
            for (unsigned k {0}; k < K; ++k) {
                result += taus[s][n][k] * log_likelihoods[s][g][k][n];
            }
        }
    }
    return result;
}

// E[ln p(R | Z, g)]
template <std::size_t K>
auto expectation(const ProbabilityVector& genotype_posteriors,
                   const VBResponsabilityMatrix<K>& taus,
                   const VBReadLikelihoodMatrix<K>& log_likelihoods)
{
    double result {0};
    for (std::size_t g {0}; g < genotype_posteriors.size(); ++g) {
        result += genotype_posteriors[g] * expectation(taus, log_likelihoods, g);
    }
    return result;
}

template <std::size_t K>
auto dirichlet_expectation(const VBAlpha<K>& posterior)
{
    using boost::math::digamma;
    const auto da0 = digamma(sum(posterior));
    return std::accumulate(std::cbegin(posterior), std::cend(posterior), 0.0,
                           [da0](const auto curr, const auto a) {
                               return curr + ((a - 1) * (digamma(a) - da0));
                           }) - maths::log_beta(posterior);
}

template <std::size_t K>
auto expectation(const VBAlphaVector<K>& posteriors)
{
    return std::accumulate(std::cbegin(posteriors), std::cend(posteriors), 0.0,
                           [](const auto curr, const auto& posterior) {
                               return curr + dirichlet_expectation(posterior);
                           });
}

template <std::size_t K>
auto q_expectation(const VBTau<K>& tau) noexcept
{
    return std::accumulate(std::cbegin(tau), std::cend(tau), 0.0,
                           [](const auto curr, const auto t) noexcept {
                               return curr + (t * std::log(t));
                           });
}

template <>
inline auto q_expectation<2>(const VBTau<2>& tau) noexcept
{
    return tau[0] * std::log(tau[0]) + tau[1] * std::log(tau[1]);
}

// E [ln q(Z_s)]
template <std::size_t K>
auto q_expectation(const VBResponsabilityVector<K>& taus) noexcept
{
    return std::accumulate(std::cbegin(taus), std::cend(taus), 0.0,
                           [](const auto curr, const auto& tau) noexcept {
                               return curr + q_expectation(tau);
                           });
}

// sum s E [ln q(Z_s)]
template <std::size_t K>
auto q_expectation(const VBResponsabilityMatrix<K>& taus) noexcept
{
    return std::accumulate(std::cbegin(taus), std::cend(taus), 0.0,
                           [](const auto curr, const auto& t) noexcept {
                               return curr + q_expectation(t);
                           });
}

template <std::size_t K>
auto calculate_lower_bound(const VBAlphaVector<K>& prior_alphas,
                             const LogProbabilityVector& genotype_log_priors,
                             const VBReadLikelihoodMatrix<K>& log_likelihoods,
                             const VBLatents<K>& latents)
{
    const auto& genotype_posteriors     = latents.genotype_posteriors;
    const auto& genotype_log_posteriors = latents.genotype_log_posteriors;
    const auto& posterior_alphas        = latents.alphas;
    const auto& taus                    = latents.responsabilities;
    double result {0};
    result += expectation(genotype_posteriors, genotype_log_priors);
    result += expectation(prior_alphas, posterior_alphas);
    result += expectation(taus, posterior_alphas);
    result += expectation(genotype_posteriors, taus, log_likelihoods);
    result -= expectation(genotype_posteriors, genotype_log_posteriors);
    result -= expectation(posterior_alphas);
    result -= q_expectation(taus);
    return result;
}

// Main algorithm - single seed

// Starting iteration with given genotype_log_posteriors
template <std::size_t K>
VBLatents<K>
run_variational_bayes(const VBAlphaVector<K>& prior_alphas,
                      const LogProbabilityVector& genotype_log_priors,
                      const VBReadLikelihoodMatrix<K>& log_likelihoods,
                      LogProbabilityVector genotype_log_posteriors,
                      const VariationalBayesParameters& params)
{
    assert(!prior_alphas.empty());
    assert(!genotype_log_priors.empty());
    assert(!log_likelihoods.empty());
    assert(prior_alphas.size() == log_likelihoods.size()); // num samples
    assert(log_likelihoods.front().size() == genotype_log_priors.size()); // num genotypes
    assert(params.max_iterations > 0);
    auto genotype_posteriors = exp(genotype_log_posteriors);
    auto posterior_alphas = prior_alphas;
    auto responsabilities = init_responsabilities<K>(posterior_alphas, genotype_posteriors, log_likelihoods);
    assert(responsabilities.size() == log_likelihoods.size()); // num samples
    bool is_converged {false};
    double max_change {0};
    for (unsigned i {0}; i < params.max_iterations; ++i) {
        update_genotype_log_posteriors(genotype_log_posteriors, genotype_log_priors, responsabilities, log_likelihoods);
        exp(genotype_log_posteriors, genotype_posteriors);
        update_alphas(posterior_alphas, prior_alphas, responsabilities);
        update_responsabilities(responsabilities, posterior_alphas, genotype_posteriors, log_likelihoods);
        std::tie(is_converged, max_change) = check_convergence(prior_alphas, posterior_alphas, max_change, params.epsilon);
        if (is_converged) break;
    }
    return VBLatents<K> {
        std::move(genotype_posteriors), std::move(genotype_log_posteriors),
        std::move(posterior_alphas), std::move(responsabilities)
    };
}

template <std::size_t K>
std::pair<VBLatents<K>, double>
get_max_evidence_seed(const VBAlphaVector<K>& prior_alphas,
                      const LogProbabilityVector& genotype_log_priors,
                      const VBReadLikelihoodMatrix<K>& log_likelihoods,
                      std::vector<VBLatents<K>>&& latents)
{
    std::vector<double> seed_evidences(latents.size());
    std::transform(std::cbegin(latents), std::cend(latents), std::begin(seed_evidences),
                   [&](const auto& seed_latents) {
                       return calculate_lower_bound(prior_alphas, genotype_log_priors, log_likelihoods, seed_latents);
                   });
    const auto max_itr = std::max_element(std::cbegin(seed_evidences), std::cend(seed_evidences));
    const auto max_idx = std::distance(std::cbegin(seed_evidences), max_itr);
    return std::make_pair(std::move(latents[max_idx]), *max_itr);
}

} // namespace detail

template <std::size_t K>
std::pair<VBLatents<K>, double>
run_variational_bayes(const VBAlphaVector<K>& prior_alphas,
                      const LogProbabilityVector& genotype_log_priors,
                      const VBReadLikelihoodMatrix<K>& log_likelihoods,
                      const VariationalBayesParameters& params,
                      std::vector<LogProbabilityVector> seeds)
{
    std::vector<VBLatents<K>> latents {};
    latents.reserve(seeds.size());
    for (auto& seed : seeds) {
        latents.push_back(detail::run_variational_bayes(prior_alphas, genotype_log_priors, log_likelihoods,
                                                        std::move(seed), params));
    }
    return detail::get_max_evidence_seed(prior_alphas, genotype_log_priors, log_likelihoods, std::move(latents));
}

inline VBReadLikelihoodArray::VBReadLikelihoodArray(const BaseType& underlying_likelihoods)
: likelihoods{std::addressof(underlying_likelihoods)} {}

inline void VBReadLikelihoodArray::operator=(const BaseType& other)
{
    likelihoods = std::addressof(other);
}

inline void VBReadLikelihoodArray::operator=(std::reference_wrapper<const BaseType> other)
{
    likelihoods = std::addressof(other.get());
}

inline std::size_t VBReadLikelihoodArray::size() const noexcept
{
    return likelihoods->size();
}

inline VBReadLikelihoodArray::BaseType::const_iterator VBReadLikelihoodArray::begin() const noexcept
{
    return likelihoods->begin();
}

inline VBReadLikelihoodArray::BaseType::const_iterator VBReadLikelihoodArray::end() const noexcept
{
    return likelihoods->end();
}

inline double VBReadLikelihoodArray::operator[](const std::size_t n) const noexcept
{
    return likelihoods->operator[](n);
}

} // namespace model
} // namespace octopus

#endif
