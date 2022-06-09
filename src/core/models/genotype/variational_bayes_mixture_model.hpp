// Copyright (c) 2015-2021 Daniel Cooke
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
#include <type_traits>

#include <boost/optional.hpp>
#include <boost/math/special_functions/digamma.hpp>

#include "core/models/haplotype_likelihood_array.hpp"
#include "utils/maths.hpp"
#include "utils/memory_footprint.hpp"
#include "utils/parallel_transform.hpp"
#include "utils/thread_pool.hpp"


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
    double epsilon = 0.05;
    unsigned max_iterations = 1000;
    bool save_memory = false;
};

using ProbabilityVector    = std::vector<double>;
using LogProbabilityVector = std::vector<double>;

template <std::size_t K>
using VBAlpha = std::array<float, K>;
template <std::size_t K>
using VBAlphaVector = std::vector<VBAlpha<K>>;

class VBReadLikelihoodArray
{
public:
    using BaseType = HaplotypeLikelihoodArray::LikelihoodVector;
    
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
    BaseType::value_type operator[](const std::size_t n) const noexcept;

private:
    const BaseType* likelihoods;
};

template <std::size_t K>
using VBGenotype = std::array<VBReadLikelihoodArray, K>; // One element per haplotype in genotype (i.e. K)
template <std::size_t K>
using VBGenotypeVector = std::vector<VBGenotype<K>>; // Per element per genotype
template <std::size_t K>
using VBReadLikelihoodMatrix = std::vector<VBGenotypeVector<K>>; // One element per sample

using VBTau = std::vector<VBReadLikelihoodArray::BaseType::value_type>; // One element per read
template <std::size_t K>
using VBResponsibilityVector = std::array<VBTau, K>; // One element per haplotype in genotype (i.e. K)
template <std::size_t K>
using VBResponsibilityMatrix = std::vector<VBResponsibilityVector<K>>; // One element per sample

template <std::size_t K>
struct VBLatents
{
    ProbabilityVector genotype_posteriors;
    LogProbabilityVector genotype_log_posteriors;
    VBAlphaVector<K> alphas;
    VBResponsibilityMatrix<K> responsibilities;
};

// Main VB method

namespace detail {

using VBExpandedLikelihood = std::vector<float>; // One element per genotype
using VBExpandedGenotype = std::vector<VBExpandedLikelihood>; // One element per read
template <std::size_t K>
using VBExpandedGenotypeVector = std::array<VBExpandedGenotype, K>; // One element per haplotype in genotype
template <std::size_t K>
using VBExpandedLikelihoodMatrix = std::vector<VBExpandedGenotypeVector<K>>; // One element per sample

template <std::size_t K>
auto invert(const VBGenotypeVector<K>& likelihoods)
{
    static_assert(K > 0, "K == 0");
    const auto num_genotypes = likelihoods.size();
    assert(num_genotypes > 0);
    const auto num_reads = likelihoods.front().front().size();
    VBExpandedGenotypeVector<K> result {};
    for (std::size_t k {0}; k < K; ++k) {
        result[k] = VBExpandedGenotype(num_reads, VBExpandedLikelihood(num_genotypes));
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
    VBExpandedLikelihoodMatrix<K> result {};
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
    using T = typename VBAlpha<K>::value_type;
    return std::accumulate(std::cbegin(alpha), std::cend(alpha), T {0});
}

template <typename T, std::size_t K>
auto compute_digamma_diffs(const std::array<T, K>& alphas)
{
    std::array<T, K> result;
    using boost::math::digamma;
    const auto digamma_a0 = digamma(sum(alphas));
    for (unsigned k {0}; k < K; ++k) {
        result[k] = digamma(alphas[k]) - digamma_a0;
    }
    return result;
}

template <std::size_t K>
auto count_reads(const VBGenotypeVector<K>& likelihoods) noexcept
{
    return likelihoods[0][0].size();
}

template <std::size_t K>
auto count_reads(const VBExpandedGenotypeVector<K>& likelihoods) noexcept
{
    return likelihoods[0].size();
}

template <typename ProbabilityVector_, std::size_t K>
auto marginalise(const ProbabilityVector_& distribution, const VBGenotypeVector<K>& likelihoods,
                 const unsigned k, const std::size_t n) noexcept
{
    using T = typename ProbabilityVector_::value_type;
    return std::inner_product(std::cbegin(distribution), std::cend(distribution),
                              std::cbegin(likelihoods), T {0}, std::plus<> {},
                              [k, n] (const auto p, const auto& haplotype_likelihoods) noexcept {
                                  return p * haplotype_likelihoods[k][n];
                              });
}

template <typename T1, typename T2>
auto inner_product(const T1& lhs, const T2& rhs) noexcept
{
    assert(std::distance(std::cbegin(lhs), std::cend(lhs)) == std::distance(std::cbegin(rhs), std::cend(rhs)));
    using T = typename T1::value_type;
    return std::inner_product(std::cbegin(lhs), std::cend(lhs), std::cbegin(rhs), T {0});
}

template <typename ProbabilityVector_, std::size_t K>
auto marginalise(const ProbabilityVector_& distribution, const VBExpandedGenotypeVector<K>& likelihoods,
                 const unsigned k, const std::size_t n) noexcept
{
    return inner_product(distribution, likelihoods[k][n]);
}

template <std::size_t K, typename T, typename ProbabilityVector_, typename VBLikelihoodGenotypeVector>
void
update_responsibilities_helper(VBResponsibilityVector<K>& result,
                               const std::array<T, K>& al,
                               const ProbabilityVector_& genotype_probabilities,
                               const VBLikelihoodGenotypeVector& read_likelihoods,
                               std::true_type)
{
    const auto N = count_reads(read_likelihoods);
    std::array<T, K> ln_rho;
    for (std::size_t n {0}; n < N; ++n) {
        for (unsigned k {0}; k < K; ++k) {
            ln_rho[k] = al[k] + marginalise(genotype_probabilities, read_likelihoods, k, n);
        }
        const auto ln_rho_norm = maths::fast_log_sum_exp(ln_rho);
        for (unsigned k {0}; k < K; ++k) {
            result[k][n] = maths::fast_exp(ln_rho[k] - ln_rho_norm);
        }
    }
}

template <std::size_t K, typename T, typename ProbabilityVector_>
void
update_responsibilities_helper(VBResponsibilityVector<K>& result,
                               const std::array<T, K>& al,
                               const ProbabilityVector_& genotype_probabilities,
                               const VBExpandedGenotypeVector<K>& read_likelihoods,
                               std::false_type)
{
    using LikelihoodType = VBExpandedLikelihood::value_type;
    const std::vector<LikelihoodType> demoted_genotype_probabilities {std::cbegin(genotype_probabilities), std::cend(genotype_probabilities)};
    update_responsibilities_helper(result, al, demoted_genotype_probabilities, read_likelihoods, std::true_type {});
}

template <std::size_t K, typename T, typename ProbabilityVector_>
void
update_responsibilities_helper(VBResponsibilityVector<K>& result,
                               const std::array<T, K>& al,
                               const ProbabilityVector_& genotype_probabilities,
                               const VBExpandedGenotypeVector<K>& read_likelihoods)
{
    // For the expanded likelihood array, the inner product between likelihoods and genotype
    // posteriors - a key bottleneck in the responsibility update calculate - can be vectorised.
    // To ensure optimal execution the floating point types of the genotype probabilities and likelihoods should match.
    using ProbabilityType = typename ProbabilityVector_::value_type;
    using LikelihoodType = VBExpandedLikelihood::value_type;
    update_responsibilities_helper(result, al, genotype_probabilities, read_likelihoods,
                                   std::is_same<ProbabilityType, LikelihoodType> {});
}

template <std::size_t K, typename T, typename VBLikelihoodGenotypeVector, typename ProbabilityVector_>
void
update_responsibilities_helper(VBResponsibilityVector<K>& result,
                               const std::array<T, K>& al,
                               const ProbabilityVector_& genotype_probabilities,
                               const VBLikelihoodGenotypeVector& read_likelihoods)
{
    update_responsibilities_helper(result, al, genotype_probabilities, read_likelihoods, std::true_type {});
}

template <std::size_t K, typename VBLikelihoodGenotypeVector>
void update_responsibilities(VBResponsibilityVector<K>& result,
                             const VBAlpha<K>& posterior_alphas,
                             const ProbabilityVector& genotype_probabilities,
                             const VBLikelihoodGenotypeVector& read_likelihoods)
{
    const auto al = compute_digamma_diffs(posterior_alphas);
    update_responsibilities_helper(result, al, genotype_probabilities, read_likelihoods);
}

template <std::size_t K, typename VBLikelihoodMatrix>
void update_responsibilities(VBResponsibilityMatrix<K>& result,
                             const VBAlphaVector<K>& posterior_alphas,
                             const ProbabilityVector& genotype_probabilities,
                             const VBLikelihoodMatrix& read_likelihoods)
{
    const auto S = read_likelihoods.size();
    for (std::size_t s {0}; s < S; ++s) {
        update_responsibilities(result[s], posterior_alphas[s], genotype_probabilities, read_likelihoods[s]);
    }
}

template <std::size_t K, typename VBLikelihoodVector_>
VBResponsibilityVector<K>
init_responsibilities(const VBAlpha<K>& prior_alphas,
                      const ProbabilityVector& genotype_probabilities,
                      const VBLikelihoodVector_& read_likelihoods)
{
    const auto N = count_reads(read_likelihoods);
    VBResponsibilityVector<K> result {};
    for (auto& tau : result) tau.resize(N);
    update_responsibilities(result, prior_alphas, genotype_probabilities, read_likelihoods);
    return result;
}

template <std::size_t K, typename VBLikelihoodMatrix>
VBResponsibilityMatrix<K>
init_responsibilities(const VBAlphaVector<K>& prior_alphas,
                      const ProbabilityVector& genotype_probabilities,
                      const VBLikelihoodMatrix& read_likelihoods)
{
    const auto S = read_likelihoods.size(); // num samples
    VBResponsibilityMatrix<K> result {};
    result.reserve(S);
    for (std::size_t s {0}; s < S; ++s) {
        result.push_back(init_responsibilities(prior_alphas[s], genotype_probabilities, read_likelihoods[s]));
    }
    return result;
}

template <typename T>
inline auto sum(const std::vector<T>& values) noexcept
{
    return std::accumulate(std::cbegin(values), std::cend(values), T {});
}

template <std::size_t K>
void update_alpha(VBAlpha<K>& alpha, const VBAlpha<K>& prior_alpha,
                  const VBResponsibilityVector<K>& taus) noexcept
{
    for (unsigned k {0}; k < K; ++k) {
        alpha[k] = prior_alpha[k] + sum(taus[k]);
    }
}

template <std::size_t K>
void update_alphas(VBAlphaVector<K>& alphas, const VBAlphaVector<K>& prior_alphas,
                   const VBResponsibilityMatrix<K>& responsibilities) noexcept
{
    const auto S = alphas.size();
    assert(S == prior_alphas.size() && S == responsibilities.size());
    for (std::size_t s {0}; s < S; ++s) {
        update_alpha(alphas[s], prior_alphas[s], responsibilities[s]);
    }
}

inline auto marginalise(const VBTau& responsibilities, const VBReadLikelihoodArray& likelihoods) noexcept
{
    assert(responsibilities.size() == likelihoods.size()); // num reads
    return inner_product(responsibilities, likelihoods);
}

template <std::size_t K>
auto marginalise(const VBResponsibilityVector<K>& responsibilities,
                 const VBGenotype<K>& read_likelihoods) noexcept
{
    double result {0};
    for (unsigned k {0}; k < K; ++k) {
        result += marginalise(responsibilities[k], read_likelihoods[k]);
    }
    return result;
}

template <std::size_t K>
auto marginalise(const VBResponsibilityMatrix<K>& responsibilities,
                 const VBReadLikelihoodMatrix<K>& read_likelihoods,
                 const std::size_t g) noexcept
{
    double result {0};
    const auto S = read_likelihoods.size(); // num samples
    assert(S == responsibilities.size());
    for (std::size_t s {0}; s < S; ++s) {
        result += marginalise(responsibilities[s], read_likelihoods[s][g]);
    }
    return result;
}

template <std::size_t K>
void update_genotype_log_posteriors(LogProbabilityVector& result,
                                    const LogProbabilityVector& genotype_log_priors,
                                    const VBResponsibilityMatrix<K>& responsibilities,
                                    const VBReadLikelihoodMatrix<K>& read_likelihoods)
{
    const auto G = result.size();
    for (std::size_t g {0}; g < G; ++g) {
        result[g] = genotype_log_priors[g] + marginalise(responsibilities, read_likelihoods, g);
    }
    maths::normalise_logs(result);
}

inline auto entropy(const VBTau& tau) noexcept
{
    using T = VBTau::value_type;
    return -std::accumulate(std::cbegin(tau), std::cend(tau), T {0},
                            [] (const auto curr, const auto t) noexcept { return curr + (t * std::log(t)); });
}

// E [ln q(Z_s)]
template <std::size_t K>
auto sum_entropies(const VBResponsibilityVector<K>& taus) noexcept
{
    using T = VBTau::value_type;
    return std::accumulate(std::cbegin(taus), std::cend(taus), T {0},
                           [] (const auto curr, const auto& tau) noexcept { return curr + entropy(tau); });
}

template <std::size_t K>
auto calculate_evidence_lower_bound(const VBAlphaVector<K>& prior_alphas,
                                    const VBAlphaVector<K>& posterior_alphas,
                                    const LogProbabilityVector& genotype_log_priors,
                                    const ProbabilityVector& genotype_posteriors,
                                    const LogProbabilityVector& genotype_log_posteriors,
                                    const VBResponsibilityMatrix<K>& taus,
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
    auto responsibilities = init_responsibilities<K>(posterior_alphas, genotype_posteriors, log_likelihoods2);
    assert(responsibilities.size() == log_likelihoods1.size()); // num samples
    auto prev_evidence = std::numeric_limits<double>::lowest();
    for (unsigned i {0}; i < params.max_iterations; ++i) {
        update_genotype_log_posteriors(genotype_log_posteriors, genotype_log_priors, responsibilities, log_likelihoods1);
        exp(genotype_log_posteriors, genotype_posteriors);
        update_alphas(posterior_alphas, prior_alphas, responsibilities);
        auto curr_evidence = calculate_evidence_lower_bound(prior_alphas, posterior_alphas, genotype_log_priors,
                                                            genotype_posteriors, genotype_log_posteriors, responsibilities,
                                                            log_likelihoods1, 1e-10);
        if (curr_evidence <= prev_evidence || (curr_evidence - prev_evidence) < params.epsilon) break;
        prev_evidence = curr_evidence;
        update_responsibilities(responsibilities, posterior_alphas, genotype_posteriors, log_likelihoods2);
    }
    return VBLatents<K> {
        std::move(genotype_posteriors), std::move(genotype_log_posteriors),
        std::move(posterior_alphas), std::move(responsibilities)
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
    return !params.save_memory;
}

template <std::size_t K>
std::vector<VBLatents<K>>
run_variational_bayes(const VBAlphaVector<K>& prior_alphas,
                      const LogProbabilityVector& genotype_log_priors,
                      const VBReadLikelihoodMatrix<K>& log_likelihoods,
                      const VariationalBayesParameters& params,
                      std::vector<LogProbabilityVector>&& seeds,
                      boost::optional<ThreadPool&> workers = boost::none)
{
    std::vector<VBLatents<K>> result {};
    result.reserve(seeds.size());
    if (run_vb_with_matrix_inversion(log_likelihoods, params, seeds)) {
        const auto inverted_log_likelihoods = invert(log_likelihoods);
        const auto func = [&] (auto&& seed) { return detail::run_variational_bayes(prior_alphas, genotype_log_priors, log_likelihoods,
                                                                                   inverted_log_likelihoods, std::move(seed), params); };
        if (workers) {
            transform(std::make_move_iterator(std::begin(seeds)), std::make_move_iterator(std::end(seeds)),
                      std::back_inserter(result), func, *workers);
        } else {
            for (auto& seed : seeds) result.push_back(func(std::move(seed)));
        }
    } else {
        const auto func = [&] (auto&& seed) { return detail::run_variational_bayes(prior_alphas, genotype_log_priors, log_likelihoods,
                                                                                   std::move(seed), params); };
        if (workers) {
            transform(std::make_move_iterator(std::begin(seeds)), std::make_move_iterator(std::end(seeds)),
                      std::back_inserter(result), func, *workers);
        } else {
            for (auto& seed : seeds) result.push_back(func(std::move(seed)));
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
                                          latents.responsibilities, log_likelihoods);
    
}

inline void check_normalisation(ProbabilityVector& probabilities) noexcept
{
    const auto mass = std::accumulate(std::cbegin(probabilities), std::cend(probabilities), 0.0);
    if (mass > 1.0) for (auto& p : probabilities) p /= mass;
}

template <std::size_t K>
void check_normalisation(VBLatents<K>& latents) noexcept
{
    check_normalisation(latents.genotype_posteriors);
}

inline std::size_t max_element_index(const ProbabilityVector& probabilities) noexcept
{
    return std::distance(std::cbegin(probabilities), std::max_element(std::cbegin(probabilities), std::cend(probabilities)));
}

template <std::size_t K>
auto
find_map_modes(const std::vector<VBLatents<K>>& latents, 
               const std::vector<double>& log_evidences)
{
    const auto num_genotypes = latents.front().genotype_posteriors.size(); 
    std::vector<std::size_t> map_genotypes(num_genotypes, latents.size());
    for (std::size_t i {0}; i < latents.size(); ++i) {
        const auto map_genotype_idx = max_element_index(latents[i].genotype_posteriors);
        if (map_genotypes[map_genotype_idx] == latents.size() || log_evidences[i] > log_evidences[map_genotypes[map_genotype_idx]]) {
            map_genotypes[map_genotype_idx] = i;
        }
    }
    std::vector<std::size_t> result {};
    result.reserve(latents.size());
    for (std::size_t g {0}; g < num_genotypes; ++g) {
        if (map_genotypes[g] < latents.size()) {
            result.push_back(map_genotypes[g]);
        }
    }
    return result;
}

template <std::size_t K>
ProbabilityVector
compute_evidence_weighted_genotype_posteriors(const std::vector<VBLatents<K>>& latents, 
                                              const std::vector<double>& log_evidences)
{
    assert(!latents.empty());
    assert(latents.size() == log_evidences.size());
    const auto modes = find_map_modes(latents, log_evidences);
    std::vector<double> mode_weights(modes.size());
    std::transform(std::cbegin(modes), std::cend(modes), std::begin(mode_weights), [&] (auto mode) { return log_evidences[mode]; });
    maths::normalise_exp(mode_weights);
    const auto num_genotypes = latents.front().genotype_posteriors.size();    
    ProbabilityVector result(num_genotypes);
    for (std::size_t i {0}; i < modes.size(); ++i) {
        const auto mode_weight = mode_weights[i];
        const auto& mode_genotype_posteriors = latents[modes[i]].genotype_posteriors;
        std::transform(std::cbegin(mode_genotype_posteriors), std::cend(mode_genotype_posteriors),
                       std::cbegin(result), std::begin(result), 
                       [mode_weight] (auto seed_posterior, auto curr_posterior) {
                           return curr_posterior + mode_weight * seed_posterior;
                       });
    }
    check_normalisation(result);
    return result;
}

template <std::size_t K>
std::vector<double>
calculate_log_evidences(const std::vector<VBLatents<K>>& latents,
                        const VBAlphaVector<K>& prior_alphas,
                        const LogProbabilityVector& genotype_log_priors,
                        const VBReadLikelihoodMatrix<K>& log_likelihoods,
                        boost::optional<ThreadPool&> workers = boost::none)
{
    std::vector<double> result(latents.size());
    transform(std::cbegin(latents), std::cend(latents), std::begin(result),
              [&] (const auto& seed_latents) { 
                  return calculate_evidence_lower_bound(prior_alphas, genotype_log_priors, log_likelihoods, seed_latents);
              }, workers);
    return result;
}

} // namespace detail

template <std::size_t K>
struct VBResultPacket
{
    VBLatents<K> map_latents;
    double max_log_evidence;
    ProbabilityVector evidence_weighted_genotype_posteriors;
};

template <std::size_t K>
VBResultPacket<K>
run_variational_bayes(const VBAlphaVector<K>& prior_alphas,
                      const LogProbabilityVector& genotype_log_priors,
                      const VBReadLikelihoodMatrix<K>& log_likelihoods,
                      const VariationalBayesParameters& params,
                      std::vector<LogProbabilityVector> seeds,
                      boost::optional<ThreadPool&> workers = boost::none)
{
    assert(!seeds.empty());
    auto latents = detail::run_variational_bayes(prior_alphas, genotype_log_priors, log_likelihoods, params, std::move(seeds), workers);
    const auto log_evidences = detail::calculate_log_evidences(latents, prior_alphas, genotype_log_priors, log_likelihoods, workers);
    auto weighted_genotype_posteriors = detail::compute_evidence_weighted_genotype_posteriors(latents, log_evidences);
    auto max_evidence_idx = detail::max_element_index(log_evidences);
    detail::check_normalisation(latents[max_evidence_idx]);
    return {std::move(latents[max_evidence_idx]), log_evidences[max_evidence_idx], std::move(weighted_genotype_posteriors)};
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

inline VBReadLikelihoodArray::BaseType::value_type VBReadLikelihoodArray::operator[](const std::size_t n) const noexcept
{
    return likelihoods->operator[](n);
}

template <std::size_t K>
MemoryFootprint
estimate_memory_requirement(const std::vector<SampleName>& samples,
                            const HaplotypeLikelihoodArray& likelihoods,
                            const std::size_t num_genotypes,
                            VariationalBayesParameters params)
{
    std::size_t bytes {};
    for (const auto& sample : samples) {
        bytes += sizeof(VBReadLikelihoodMatrix<K>);
        bytes += sizeof(VBGenotypeVector<K>) * num_genotypes;
        bytes += sizeof(VBResponsibilityMatrix<K>);
        const auto num_likelihoods = likelihoods.num_likelihoods(sample);
        const auto tau_bytes = num_likelihoods * sizeof(VBTau::value_type);
        bytes += tau_bytes * K + sizeof(VBResponsibilityVector<K>);
        if (!params.save_memory) {
            bytes += sizeof(detail::VBExpandedLikelihoodMatrix<K>);
            auto inverse_bytes = sizeof(detail::VBExpandedLikelihood::value_type) * num_genotypes + sizeof(detail::VBExpandedLikelihood);
            inverse_bytes *= num_likelihoods;
            inverse_bytes += sizeof(detail::VBExpandedGenotype);
            bytes += K * inverse_bytes + sizeof(detail::VBExpandedGenotypeVector<K>);
        }
    }
    return MemoryFootprint {bytes};
};

} // namespace model
} // namespace octopus

#endif
