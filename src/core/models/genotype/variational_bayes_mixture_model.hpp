// Copyright (c) 2015-2018 Daniel Cooke
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
#include <limits>

#include <boost/optional.hpp>
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
using VBGenotype = std::array<VBReadLikelihoodArray, K>; // One element per haplotype in genotype (i.e. K)
template <std::size_t K>
using VBGenotypeVector = std::vector<VBGenotype<K>>; // Per element per genotype
template <std::size_t K>
using VBReadLikelihoodMatrix = std::vector<VBGenotypeVector<K>>; // One element per sample

using VBTau = std::vector<double>; // One element per read
template <std::size_t K>
using VBResponsabilityVector = std::array<VBTau, K>; // One element per haplotype in genotype (i.e. K)
template <std::size_t K>
using VBResponsabilityMatrix = std::vector<VBResponsabilityVector<K>>; // One element per sample

template <std::size_t K>
struct VBLatents
{
    ProbabilityVector genotype_posteriors;
    LogProbabilityVector genotype_log_posteriors;
    VBAlphaVector<K> alphas;
    VBResponsabilityMatrix<K> responsabilities;
};

// Main VB method

template <std::size_t K>
std::pair<VBLatents<K>, double>
run_variational_bayes(const VBAlphaVector<K>& prior_alphas,
                      const LogProbabilityVector& genotype_log_priors,
                      const VBReadLikelihoodMatrix<K>& log_likelihoods,
                      const VariationalBayesParameters& params,
                      std::vector<LogProbabilityVector> seeds);

namespace detail {

using VBInverseLikelihood = std::vector<double>; // One element per genotype
using VBInverseGenotype = std::vector<VBInverseLikelihood>; // One element per read
template <std::size_t K>
using VBInverseGenotypeVector = std::array<VBInverseGenotype, K>; // One element per haplotype in genotype
template <std::size_t K>
using VBInverseReadLikelihoodMatrix = std::vector<VBInverseGenotypeVector<K>>; // One element per sample

template <std::size_t K>
auto invert(const VBGenotypeVector<K>& likelihoods)
{
    static_assert(K > 0, "K == 0");
    const auto num_genotypes = likelihoods.size();
    assert(num_genotypes > 0);
    const auto num_reads = likelihoods.front().front().size();
    VBInverseGenotypeVector<K> result {};
    for (std::size_t k {0}; k < K; ++k) {
        result[k] = VBInverseGenotype(num_reads, VBInverseLikelihood(num_genotypes));
        for (std::size_t n {0}; n < num_reads; ++n) {
            for (std::size_t g {0}; g < num_genotypes; ++g) {
                result[k][n][g] = likelihoods[g][k][n];
            }
        }
    }
    return result;
}

template <std::size_t K>
auto invert(const VBReadLikelihoodMatrix<K>& matrix)
{
    VBInverseReadLikelihoodMatrix<K> result {};
    result.reserve(matrix.size());
    std::transform(std::cbegin(matrix), std::cend(matrix), std::back_inserter(result),
                   [] (const auto& v) { return invert(v); });
    return result;
}

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

template <typename T>
inline auto digamma_diff(const T a, const T b)
{
    using boost::math::digamma;
    return digamma(a) - digamma(b);
}

template <typename T>
inline T log_sum_exp(const std::array<T, 1>& logs)
{
    return logs[0];
}

template <typename T>
inline T log_sum_exp(const std::array<T, 2>& logs)
{
    return maths::log_sum_exp(logs[0], logs[1]);
}

template <typename T>
inline T log_sum_exp(const std::array<T, 3>& logs)
{
    return maths::log_sum_exp(logs[0], logs[1], logs[2]);
}

template <typename T, std::size_t K>
T log_sum_exp(const std::array<T, K>& logs)
{
    return maths::log_sum_exp(logs);
}

template <std::size_t K>
auto count_reads(const VBGenotypeVector<K>& likelihoods) noexcept
{
    return likelihoods[0][0].size();
}

template <std::size_t K>
auto count_reads(const VBInverseGenotypeVector<K>& likelihoods) noexcept
{
    return likelihoods[0].size();
}

template <std::size_t K>
auto marginalise(const ProbabilityVector& distribution, const VBGenotypeVector<K>& likelihoods,
                 const unsigned k, const std::size_t n) noexcept
{
    return std::inner_product(std::cbegin(distribution), std::cend(distribution),
                              std::cbegin(likelihoods), 0.0, std::plus<> {},
                              [k, n] (const auto p, const auto& haplotype_likelihoods) noexcept {
                                  return p * haplotype_likelihoods[k][n];
                              });
}

template <std::size_t K>
auto marginalise(const ProbabilityVector& distribution, const VBInverseGenotypeVector<K>& likelihoods,
                 const unsigned k, const std::size_t n) noexcept
{
    return std::inner_product(std::cbegin(distribution), std::cend(distribution), std::cbegin(likelihoods[k][n]), 0.0);
}

template <std::size_t K, typename VBLikelihoodVector_>
VBResponsabilityVector<K>
init_responsabilities(const VBAlpha<K>& prior_alphas,
                      const ProbabilityVector& genotype_probabilities,
                      const VBLikelihoodVector_& read_likelihoods)
{
    using T = typename VBAlpha<K>::value_type;
    std::array<T, K> al; // no need to keep recomputing this
    const auto a0 = sum(prior_alphas);
    for (unsigned k {0}; k < K; ++k) {
        al[k] = digamma_diff(prior_alphas[k], a0);
    }
    const auto N = count_reads(read_likelihoods);
    VBResponsabilityVector<K> result {};
    for (auto& tau : result) tau.resize(N);
    std::array<T, K> ln_rho;
    for (std::size_t n {0}; n < N; ++n) {
        for (unsigned k {0}; k < K; ++k) {
            ln_rho[k] = al[k] + marginalise(genotype_probabilities, read_likelihoods, k, n);
        }
        const auto ln_rho_norm = log_sum_exp(ln_rho);
        for (unsigned k {0}; k < K; ++k) {
            result[k][n] = std::exp(ln_rho[k] - ln_rho_norm);
        }
    }
    return result;
}

template <std::size_t K, typename VBLikelihoodMatrix>
VBResponsabilityMatrix<K>
init_responsabilities(const VBAlphaVector<K>& prior_alphas,
                      const ProbabilityVector& genotype_probabilities,
                      const VBLikelihoodMatrix& read_likelihoods)
{
    const auto S = read_likelihoods.size(); // num samples
    VBResponsabilityMatrix<K> result{};
    result.reserve(S);
    for (std::size_t s {0}; s < S; ++s) {
        result.push_back(init_responsabilities(prior_alphas[s], genotype_probabilities, read_likelihoods[s]));
    }
    return result;
}

template <std::size_t K, typename VBLikelihoodVector_>
void update_responsabilities(VBResponsabilityVector<K>& result,
                             const VBAlpha<K>& posterior_alphas,
                             const ProbabilityVector& genotype_probabilities,
                             const VBLikelihoodVector_& read_likelihoods)
{
    using T = typename VBAlpha<K>::value_type;
    std::array<T, K> al;
    const auto a0 = sum(posterior_alphas);
    for (unsigned k {0}; k < K; ++k) {
        al[k] = digamma_diff(posterior_alphas[k], a0);
    }
    const auto N = count_reads(read_likelihoods);
    std::array<T, K> ln_rho;
    for (std::size_t n {0}; n < N; ++n) {
        for (unsigned k {0}; k < K; ++k) {
            ln_rho[k] = al[k] + marginalise(genotype_probabilities, read_likelihoods, k, n);
        }
        const auto ln_rho_norm = log_sum_exp(ln_rho);
        for (unsigned k {0}; k < K; ++k) {
            result[k][n] = std::exp(ln_rho[k] - ln_rho_norm);
        }
    }
}

// same as init_responsabilities but in-place
template <std::size_t K, typename VBLikelihoodMatrix>
void update_responsabilities(VBResponsabilityMatrix<K>& result,
                             const VBAlphaVector<K>& posterior_alphas,
                             const ProbabilityVector& genotype_probabilities,
                             const VBLikelihoodMatrix& read_likelihoods)
{
    const auto S = read_likelihoods.size();
    for (std::size_t s {0}; s < S; ++s) {
        update_responsabilities(result[s], posterior_alphas[s], genotype_probabilities, read_likelihoods[s]);
    }
}

template <typename T>
inline auto sum(const std::vector<T>& values) noexcept
{
    return std::accumulate(std::cbegin(values), std::cend(values), T {});
}

template <std::size_t K>
void update_alpha(VBAlpha<K>& alpha, const VBAlpha<K>& prior_alpha,
                  const VBResponsabilityVector<K>& taus) noexcept
{
    for (unsigned k {0}; k < K; ++k) {
        alpha[k] = prior_alpha[k] + sum(taus[k]);
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

inline auto marginalise(const VBTau& responsabilities, const VBReadLikelihoodArray& likelihoods) noexcept
{
    assert(responsabilities.size() == likelihoods.size()); // num reads
    return std::inner_product(std::cbegin(responsabilities), std::cend(responsabilities), std::cbegin(likelihoods), 0.0);
}

template <std::size_t K>
auto marginalise(const VBResponsabilityVector<K>& responsabilities,
                 const VBGenotype<K>& read_likelihoods) noexcept
{
    double result {0};
    for (unsigned k {0}; k < K; ++k) {
        result += marginalise(responsabilities[k], read_likelihoods[k]);
    }
    return result;
}

template <std::size_t K>
auto marginalise(const VBResponsabilityMatrix<K>& responsabilities,
                 const VBReadLikelihoodMatrix<K>& read_likelihoods,
                 const std::size_t g) noexcept
{
    double result {0};
    const auto S = read_likelihoods.size(); // num samples
    assert(S == responsabilities.size());
    for (std::size_t s {0}; s < S; ++s) {
        result += marginalise(responsabilities[s], read_likelihoods[s][g]);
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

inline auto entropy(const VBTau& tau) noexcept
{
    return -std::accumulate(std::cbegin(tau), std::cend(tau), 0.0,
                            [] (const auto curr, const auto t) noexcept { return curr + (t * std::log(t)); });
}

// E [ln q(Z_s)]
template <std::size_t K>
auto sum_entropies(const VBResponsabilityVector<K>& taus) noexcept
{
    return std::accumulate(std::cbegin(taus), std::cend(taus), 0.0,
                           [] (const auto curr, const auto& tau) noexcept { return curr + entropy(tau); });
}

template <std::size_t K>
auto calculate_evidence_lower_bound(const VBAlphaVector<K>& prior_alphas,
                                    const VBAlphaVector<K>& posterior_alphas,
                                    const LogProbabilityVector& genotype_log_priors,
                                    const ProbabilityVector& genotype_posteriors,
                                    const LogProbabilityVector& genotype_log_posteriors,
                                    const VBResponsabilityMatrix<K>& taus,
                                    const VBReadLikelihoodMatrix<K>& log_likelihoods,
                                    const boost::optional<double> max_posterior_skip = boost::none)
{
    const auto G = genotype_log_priors.size();
    const auto S = log_likelihoods.size();
    double result {0};
    for (std::size_t g {0}; g < G; ++g) {
        if (!max_posterior_skip || genotype_posteriors[g] >= *max_posterior_skip) {
            auto w = genotype_log_priors[g] - genotype_log_posteriors[g];
            for (std::size_t s {0}; s < S; ++s) {
                w += marginalise(taus[s], log_likelihoods[s][g]);
            }
            result += genotype_posteriors[g] * w;
        }
    }
    for (std::size_t s {0}; s < S; ++s) {
        result += (maths::log_beta(posterior_alphas[s]) - maths::log_beta(prior_alphas[s]));
        result += sum_entropies(taus[s]);
    }
    return result;
}

// Main algorithm - single seed

// Starting iteration with given genotype_log_posteriors
template <std::size_t K, typename VBLikelihoodMatrix1, typename VBLikelihoodMatrix2>
VBLatents<K>
run_variational_bayes(const VBAlphaVector<K>& prior_alphas,
                      const LogProbabilityVector& genotype_log_priors,
                      const VBLikelihoodMatrix1& log_likelihoods1,
                      const VBLikelihoodMatrix2& log_likelihoods2,
                      LogProbabilityVector genotype_log_posteriors,
                      const VariationalBayesParameters& params)
{
    assert(!prior_alphas.empty());
    assert(!genotype_log_priors.empty());
    assert(!log_likelihoods1.empty());
    assert(log_likelihoods1.size() == log_likelihoods2.size());
    assert(prior_alphas.size() == log_likelihoods1.size()); // num samples
    assert(log_likelihoods1.front().size() == genotype_log_priors.size()); // num genotypes
    assert(params.max_iterations > 0);
    auto genotype_posteriors = exp(genotype_log_posteriors);
    auto posterior_alphas = prior_alphas;
    auto responsabilities = init_responsabilities<K>(posterior_alphas, genotype_posteriors, log_likelihoods2);
    assert(responsabilities.size() == log_likelihoods1.size()); // num samples
    auto prev_evidence = std::numeric_limits<double>::lowest();
    bool is_converged {};
    double max_change {};
    for (unsigned i {0}; i < params.max_iterations; ++i) {
        update_genotype_log_posteriors(genotype_log_posteriors, genotype_log_priors, responsabilities, log_likelihoods1);
        exp(genotype_log_posteriors, genotype_posteriors);
        update_alphas(posterior_alphas, prior_alphas, responsabilities);
        std::tie(is_converged, max_change) = check_convergence(prior_alphas, posterior_alphas, max_change, params.epsilon);
        if (is_converged) break;
        auto curr_evidence = calculate_evidence_lower_bound(prior_alphas, posterior_alphas, genotype_log_priors,
                                                            genotype_posteriors, genotype_log_posteriors, responsabilities,
                                                            log_likelihoods1, 1e-10);
        if (curr_evidence <= prev_evidence || (curr_evidence - prev_evidence) < params.epsilon) break;
        prev_evidence = curr_evidence;
        update_responsabilities(responsabilities, posterior_alphas, genotype_posteriors, log_likelihoods2);
    }
    return VBLatents<K> {
        std::move(genotype_posteriors), std::move(genotype_log_posteriors),
        std::move(posterior_alphas), std::move(responsabilities)
    };
}

// Not using inverted log likelihoods
template <std::size_t K>
VBLatents<K>
run_variational_bayes(const VBAlphaVector<K>& prior_alphas,
                      const LogProbabilityVector& genotype_log_priors,
                      const VBReadLikelihoodMatrix<K>& log_likelihoods,
                      LogProbabilityVector genotype_log_posteriors,
                      const VariationalBayesParameters& params)
{
    return run_variational_bayes(prior_alphas, genotype_log_posteriors, log_likelihoods,
                                 log_likelihoods, genotype_log_posteriors, params);
}

// Main algorithm - multiple seed

template <std::size_t K>
bool run_vb_with_matrix_inversion(const VBReadLikelihoodMatrix<K>& log_likelihoods,
                                  const VariationalBayesParameters& params,
                                  const std::vector<LogProbabilityVector>& seeds) noexcept
{
    return true;
}

template <std::size_t K>
std::vector<VBLatents<K>>
run_variational_bayes(const VBAlphaVector<K>& prior_alphas,
                      const LogProbabilityVector& genotype_log_priors,
                      const VBReadLikelihoodMatrix<K>& log_likelihoods,
                      const VariationalBayesParameters& params,
                      std::vector<LogProbabilityVector>&& seeds)
{
    std::vector<VBLatents<K>> result {};
    result.reserve(seeds.size());
    if (run_vb_with_matrix_inversion(log_likelihoods, params, seeds)) {
        const auto inverted_log_likelihoods = invert(log_likelihoods);
        for (auto& seed : seeds) {
            result.push_back(detail::run_variational_bayes(prior_alphas, genotype_log_priors,
                                                           log_likelihoods, inverted_log_likelihoods,
                                                           std::move(seed), params));
        }
    } else {
        for (auto& seed : seeds) {
            result.push_back(detail::run_variational_bayes(prior_alphas, genotype_log_priors,
                                                           log_likelihoods,
                                                           std::move(seed), params));
        }
    }
    return result;
}

// lower-bound calculation

template <std::size_t K>
auto calculate_evidence_lower_bound(const VBAlphaVector<K>& prior_alphas,
                                    const LogProbabilityVector& genotype_log_priors,
                                    const VBReadLikelihoodMatrix<K>& log_likelihoods,
                                    const VBLatents<K>& latents)
{
    return calculate_evidence_lower_bound(prior_alphas, latents.alphas, genotype_log_priors,
                                          latents.genotype_posteriors, latents.genotype_log_posteriors,
                                          latents.responsabilities, log_likelihoods);
    
}

template <std::size_t K>
std::pair<VBLatents<K>, double>
get_max_evidence_latents(const VBAlphaVector<K>& prior_alphas,
                         const LogProbabilityVector& genotype_log_priors,
                         const VBReadLikelihoodMatrix<K>& log_likelihoods,
                         std::vector<VBLatents<K>>&& latents)
{
    std::vector<double> seed_evidences(latents.size());
    std::transform(std::cbegin(latents), std::cend(latents), std::begin(seed_evidences),
                   [&] (const auto& seed_latents) {
                       return calculate_evidence_lower_bound(prior_alphas, genotype_log_priors, log_likelihoods, seed_latents);
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
    assert(!seeds.empty());
    auto latents = detail::run_variational_bayes(prior_alphas, genotype_log_priors, log_likelihoods, params, std::move(seeds));
    return detail::get_max_evidence_latents(prior_alphas, genotype_log_priors, log_likelihoods, std::move(latents));
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
