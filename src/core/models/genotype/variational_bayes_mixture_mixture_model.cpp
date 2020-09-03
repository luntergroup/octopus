// Copyright (c) 2015-2020 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "variational_bayes_mixture_mixture_model.hpp"

#include <utility>
#include <cmath>
#include <iterator>
#include <algorithm>
#include <numeric>
#include <iostream>
#include <limits>

#include <boost/math/special_functions/digamma.hpp>

#include "utils/maths.hpp"

namespace octopus { namespace model {

VariationalBayesMixtureMixtureModel::VariationalBayesMixtureMixtureModel(Options options)
: options_ {std::move(options)}
{}

VariationalBayesMixtureMixtureModel::Inferences
VariationalBayesMixtureMixtureModel::evaluate(const LogProbabilityVector& genotype_log_priors,
                                              const HaplotypeLikelihoodMatrix& log_likelihoods,
                                              const GroupOptionalPriorArray& group_priors,
                                              const GroupConcentrationVector& group_concentrations,
                                              const MixtureConcentrationArray& mixture_concentrations,
                                              std::vector<LogProbabilityVector> seeds) const
{
    const auto group_log_priors = to_logs(group_priors);
    const auto evaluate_seed = [&] (auto&& seed) {
        return this->evaluate(genotype_log_priors, log_likelihoods, group_log_priors, group_concentrations, mixture_concentrations, std::move(seed)); };
    std::vector<PointInferences> seed_inferences(seeds.size());
    if (options_.parallel_execution) {
        parallel_transform(std::make_move_iterator(std::begin(seeds)), std::make_move_iterator(std::end(seeds)), std::begin(seed_inferences), evaluate_seed);
    } else {
        std::transform(std::make_move_iterator(std::begin(seeds)), std::make_move_iterator(std::end(seeds)), std::begin(seed_inferences), evaluate_seed);
    }
    Inferences result {};
    compute_evidence_weighted_latents(result.weighted_genotype_posteriors, result.weighted_group_responsibilities, seed_inferences);
    const static auto evidence_less = [] (const auto& lhs, const auto& rhs) { return lhs.approx_log_evidence < rhs.approx_log_evidence; };
    result.map = *std::max_element(std::begin(seed_inferences), std::end(seed_inferences), evidence_less);
    return result;
}

VariationalBayesMixtureMixtureModel::Inferences
VariationalBayesMixtureMixtureModel::evaluate(const LogProbabilityVector& genotype_log_priors,
                                              const HaplotypeLikelihoodMatrix& log_likelihoods,
                                              const GroupConcentrationVector& group_concentrations,
                                              const MixtureConcentrationArray& mixture_concentrations,
                                              std::vector<LogProbabilityVector> seeds) const
{
    const auto num_samples = log_likelihoods.size();
    const GroupOptionalPriorArray no_priors(num_samples);
    return evaluate(genotype_log_priors, log_likelihoods, no_priors, group_concentrations, mixture_concentrations, std::move(seeds));
}

VariationalBayesMixtureMixtureModel::Inferences
VariationalBayesMixtureMixtureModel::evaluate(const LogProbabilityVector& genotype_log_priors,
                                              const HaplotypeLikelihoodMatrix& log_likelihoods,
                                              const GroupOptionalPriorArray& group_priors,
                                              const double group_concentration,
                                              const std::vector<double>& mixture_concentrations,
                                              std::vector<LogProbabilityVector> seeds) const
{
    const auto num_samples = log_likelihoods.size();
    assert(mixture_concentrations.size() == num_samples);
    const auto num_groups = log_likelihoods.front().front().size();
    std::vector<unsigned> group_mixture_sizes {};
    group_mixture_sizes.reserve(num_groups);
    for (const auto& genotype : log_likelihoods.front().front()) {
        group_mixture_sizes.push_back(genotype.size());
    }
    const GroupConcentrationVector group_concentrations(num_groups, group_concentration);
    MixtureConcentrationArray group_mixture_concentrations {};
    group_mixture_concentrations.reserve(num_samples);
    for (auto sample_concentration : mixture_concentrations) {
        group_mixture_concentrations.emplace_back(num_groups, ComponentConcentrationVector(sample_concentration));
        unsigned t {0};
        for (auto& sample_group_sample_concentrations : group_mixture_concentrations.back()) {
            sample_group_sample_concentrations.assign(group_mixture_sizes[t++], sample_concentration);
        }
    }
    return evaluate(genotype_log_priors, log_likelihoods, group_priors, group_concentrations, group_mixture_concentrations, std::move(seeds));
}

VariationalBayesMixtureMixtureModel::Inferences
VariationalBayesMixtureMixtureModel::evaluate(const LogProbabilityVector& genotype_log_priors,
                                              const HaplotypeLikelihoodMatrix& log_likelihoods,
                                              const double group_concentration,
                                              const std::vector<double>& mixture_concentrations,
                                              std::vector<LogProbabilityVector> seeds) const
{
    const auto num_samples = log_likelihoods.size();
    const GroupOptionalPriorArray no_priors(num_samples);
    return evaluate(genotype_log_priors, log_likelihoods, no_priors, group_concentration, mixture_concentrations, std::move(seeds));
}

VariationalBayesMixtureMixtureModel::Inferences
VariationalBayesMixtureMixtureModel::evaluate(const LogProbabilityVector& genotype_log_priors,
                                              const HaplotypeLikelihoodMatrix& log_likelihoods,
                                              const GroupOptionalPriorArray& group_priors,
                                              const double group_concentration,
                                              const double mixture_concentration,
                                              std::vector<LogProbabilityVector> seeds) const
{
    const auto num_samples = log_likelihoods.size();
    const auto num_groups = log_likelihoods.front().front().size();
    std::vector<unsigned> group_mixture_sizes {};
    group_mixture_sizes.reserve(num_groups);
    for (const auto& genotype : log_likelihoods.front().front()) {
        group_mixture_sizes.push_back(genotype.size());
    }
    const GroupConcentrationVector group_concentrations(num_groups, group_concentration);
    MixtureConcentrationArray mixture_concentrations(num_samples, MixtureConcentrationVector(num_groups));
    for (auto& sample_concentrations : mixture_concentrations) {
        unsigned t {0};
        for (auto& sample_group_sample_concentrations : sample_concentrations) {
            sample_group_sample_concentrations.assign(group_mixture_sizes[t++], mixture_concentration);
        }
    }
    return evaluate(genotype_log_priors, log_likelihoods, group_priors, group_concentrations, mixture_concentrations, std::move(seeds));
}

VariationalBayesMixtureMixtureModel::Inferences
VariationalBayesMixtureMixtureModel::evaluate(const LogProbabilityVector& genotype_log_priors,
                                              const HaplotypeLikelihoodMatrix& log_likelihoods,
                                              const double group_concentration,
                                              const double mixture_concentration,
                                              std::vector<LogProbabilityVector> seeds) const
{
    const auto num_samples = log_likelihoods.size();
    const GroupOptionalPriorArray no_priors(num_samples);
    return evaluate(genotype_log_priors, log_likelihoods, no_priors, group_concentration, mixture_concentration, std::move(seeds));
}

// Private methods

VariationalBayesMixtureMixtureModel::GroupOptionalLogPriorVector
VariationalBayesMixtureMixtureModel::to_logs(const GroupOptionalPriorVector& prior) const
{
    GroupOptionalLogPriorVector result {};
    if (prior) {
        result = LogProbabilityVector(prior->size());
        const static auto to_log = [] (auto p) noexcept {
            const static auto min_log_prior = std::log(std::numeric_limits<LogProbability>::min());
            return p > 0 ? std::log(p) : min_log_prior;
        };
        std::transform(std::cbegin(*prior), std::cend(*prior), std::begin(*result), to_log);
    }
    return result;
}

VariationalBayesMixtureMixtureModel::GroupOptionalLogPriorArray
VariationalBayesMixtureMixtureModel::to_logs(const GroupOptionalPriorArray& priors) const
{
    GroupOptionalLogPriorArray result {};
    result.reserve(priors.size());
    std::transform(std::cbegin(priors), std::cend(priors), std::back_inserter(result),
                   [this] (const auto& sample_priors) { return to_logs(sample_priors); });
    return result;
}

namespace {

VariationalBayesMixtureMixtureModel::ProbabilityVector&
exp(const LogProbabilityVector& log_probabilities, ProbabilityVector& result) noexcept
{
    std::transform(std::cbegin(log_probabilities), std::cend(log_probabilities), std::begin(result),
                   [] (const auto lp) noexcept { return std::exp(lp); });
    return result;
}

auto exp(const VariationalBayesMixtureMixtureModel::LogProbabilityVector& log_probabilities)
{
    VariationalBayesMixtureMixtureModel::ProbabilityVector result(log_probabilities.size());
    return exp(log_probabilities, result);
}

template <typename T>
auto sum(const std::vector<T>& values) noexcept
{
    return std::accumulate(std::cbegin(values), std::cend(values), T {0});
}

auto sum(const VBReadLikelihoodArray& likelihoods) noexcept
{
    using T = VBReadLikelihoodArray::BaseType::value_type;
    return std::accumulate(std::cbegin(likelihoods), std::cend(likelihoods), T {0});
}

template <typename T1, typename T2>
auto inner_product(const T1& lhs, const T2& rhs) noexcept
{
    assert(std::distance(std::cbegin(lhs), std::cend(lhs)) == std::distance(std::cbegin(rhs), std::cend(rhs)));
    using T = typename T1::value_type;
    return std::inner_product(std::cbegin(lhs), std::cend(lhs), std::cbegin(rhs), T {0});
}

template <typename T>
inline auto digamma_diff(const T a, const T b)
{
    using boost::math::digamma;
    return digamma(a) - digamma(b);
}

template <typename T>
auto dirichlet_expectation_log(const std::vector<T>& concentrations)
{
    std::vector<T> result(concentrations.size());
    if (concentrations.size() > 1) {
        const auto alpha_0 = sum(concentrations);
        std::transform(std::cbegin(concentrations), std::cend(concentrations), std::begin(result),
                       [=] (auto alpha) { return digamma_diff(alpha, alpha_0); });
    }
    return result;
}

template <typename T>
T log_sum_exp(const std::vector<T>& logs)
{
    return maths::log_sum_exp(logs);
}

bool all_equal_sizes(const VariationalBayesMixtureMixtureModel::MixtureConcentrationVector& concentrations) noexcept
{
    const auto size_unequal = [] (const auto& lhs, const auto& rhs) { return lhs.size() != rhs.size(); };
    return std::adjacent_find(std::cbegin(concentrations), std::end(concentrations), size_unequal) == std::cend(concentrations);
}

} // namespace

VariationalBayesMixtureMixtureModel::PointInferences
VariationalBayesMixtureMixtureModel::evaluate(const LogProbabilityVector& genotype_log_priors,
                                              const HaplotypeLikelihoodMatrix& log_likelihoods,
                                              const GroupOptionalLogPriorArray& group_log_priors,
                                              const GroupConcentrationVector& prior_group_concentrations,
                                              const MixtureConcentrationArray& prior_mixture_concentrations,
                                              LogProbabilityVector genotype_log_posteriors) const
{
    Latents latents {};
    latents.genotype_log_posteriors = std::move(genotype_log_posteriors);
    latents.genotype_posteriors = exp(latents.genotype_log_posteriors);
    latents.group_concentrations = prior_group_concentrations;
    latents.mixture_concentrations = prior_mixture_concentrations;
    latents.group_responsibilities = init_responsibilities(group_log_priors, prior_group_concentrations, prior_mixture_concentrations,
                                                           latents.genotype_posteriors, log_likelihoods);
    latents.component_responsibilities = init_responsibilities(prior_group_concentrations, prior_mixture_concentrations,
                                                               latents.genotype_posteriors, latents.group_responsibilities, log_likelihoods);
    auto prev_evidence = std::numeric_limits<double>::lowest();
    for (unsigned i {0}; i < options_.max_iterations; ++i) {
        update_genotype_log_posteriors(latents.genotype_log_posteriors, genotype_log_priors,
                                       latents.group_responsibilities, latents.component_responsibilities,
                                       log_likelihoods);
        exp(latents.genotype_log_posteriors, latents.genotype_posteriors);
        update_group_concentrations(latents.group_concentrations, prior_group_concentrations, latents.group_responsibilities);
        update_mixture_concentrations(latents.mixture_concentrations, prior_mixture_concentrations,
                                      latents.group_responsibilities, latents.component_responsibilities);
        auto curr_evidence = calculate_evidence(prior_group_concentrations, latents.group_concentrations,
                                                prior_mixture_concentrations, latents.mixture_concentrations,
                                                genotype_log_priors, latents.genotype_log_posteriors, latents.genotype_posteriors,
                                                latents.group_responsibilities, latents.component_responsibilities,
                                                log_likelihoods);
//        assert(curr_evidence + options_.epsilon >= prev_evidence);
        if (curr_evidence <= prev_evidence || (curr_evidence - prev_evidence) < options_.epsilon) {
            prev_evidence = curr_evidence;
            break;
        }
        prev_evidence = curr_evidence;
        update_responsibilities(latents.group_responsibilities, group_log_priors,latents.group_concentrations,
                                latents.mixture_concentrations, latents.genotype_posteriors,
                                latents.component_responsibilities, log_likelihoods);
        update_responsibilities(latents.component_responsibilities, latents.group_concentrations, latents.mixture_concentrations,
                                latents.genotype_posteriors, latents.group_responsibilities, log_likelihoods);
    }
    return {std::move(latents), prev_evidence};
}

namespace {

template <typename ForwardIterator>
std::size_t max_element_index(ForwardIterator first, ForwardIterator last)
{
    return std::distance(first, std::max_element(first, last));
}

template <typename Range>
std::size_t max_element_index(const Range& values)
{
    return max_element_index(std::cbegin(values), std::cend(values));
}

} // namespace

VariationalBayesMixtureMixtureModel::GroupResponsibilityVector
VariationalBayesMixtureMixtureModel::init_responsibilities(const GroupOptionalLogPriorArray& group_log_priors,
                                                           const GroupConcentrationVector& group_concentrations,
                                                           const MixtureConcentrationArray& mixture_concentrations,
                                                           const ProbabilityVector& genotype_priors,
                                                           const HaplotypeLikelihoodMatrix& log_likelihoods) const
{
    const auto S = log_likelihoods.size();
    const auto T = group_concentrations.size();
    GroupResponsibilityVector result(S, Sigma(T));
    const auto G = genotype_priors.size();
    const auto ln_ex_psi = dirichlet_expectation_log(group_concentrations);
    std::size_t max_K {0};
    for (std::size_t s {0}; s < S; ++s) {
        for (std::size_t t {0}; t < T; ++t) {
            max_K = std::max(max_K, mixture_concentrations[s][t].size());
        }
    }
    const auto max_genotype_prior_idx = max_element_index(genotype_priors);
    ComponentResponsibilityVector approx_taus(max_K);
    for (std::size_t s {0}; s < S; ++s) {
        const auto N = log_likelihoods[s][0][0][0].size(); // num reads
        for (auto& tau : approx_taus) tau.resize(N);
        for (std::size_t t {0}; t < T; ++t) {
            result[s][t] = ln_ex_psi[t];
            if (group_log_priors[s]) result[s][t] += (*group_log_priors[s])[t];
            const auto ln_ex_pi = dirichlet_expectation_log(mixture_concentrations[s][t]);
            const auto K = mixture_concentrations[s][t].size();
            // Approximate tau assuming p(t) = 1
            Tau approx_tau(K);
            for (std::size_t n {0}; n < N; ++n) {
                for (std::size_t k {0}; k < K; ++k) {
                    approx_tau[k] = ln_ex_pi[k] + log_likelihoods[s][max_genotype_prior_idx][t][k][n];
                }
                maths::normalise_exp(approx_tau);
                for (std::size_t k {0}; k < K; ++k) {
                    approx_taus[k][n] = approx_tau[k];
                }
            }
            for (std::size_t k {0}; k < K; ++k) {
                for (std::size_t g {0}; g < G; ++g) {
                    result[s][t] += genotype_priors[g] * inner_product(approx_taus[k], log_likelihoods[s][g][t][k]);
                }
            }
        }
        maths::normalise_exp(result[s]);
    }
    return result;
}

void
VariationalBayesMixtureMixtureModel::update_responsibilities(GroupResponsibilityVector& result,
                                                             const GroupOptionalLogPriorArray& group_log_priors,
                                                             const GroupConcentrationVector& group_concentrations,
                                                             const MixtureConcentrationArray& mixture_concentrations,
                                                             const ProbabilityVector& genotype_posteriors,
                                                             const ComponentResponsibilityMatrix& component_responsibilities,
                                                             const HaplotypeLikelihoodMatrix& log_likelihoods) const
{
    const auto T = group_concentrations.size();
    const auto S = log_likelihoods.size();
    const auto G = genotype_posteriors.size();
    const auto ln_ex_psi = dirichlet_expectation_log(group_concentrations);
    for (std::size_t s {0}; s < S; ++s) {
        for (std::size_t t {0}; t < T; ++t) {
            const auto max_K = component_responsibilities[s][t].size();
            result[s][t] = ln_ex_psi[t];
            if (group_log_priors[s]) result[s][t] += (*group_log_priors[s])[t];
            const auto ln_ex_pi = dirichlet_expectation_log(mixture_concentrations[s][t]);
            const auto K = ln_ex_pi.size();
            for (std::size_t k {0}; k < max_K; ++k) {
                for (std::size_t g {0}; g < G; ++g) {
                    if (k < K) {
                        result[s][t] += genotype_posteriors[g] * inner_product(component_responsibilities[s][t][k], log_likelihoods[s][g][t][k]);
                    } else {
                        result[s][t] += genotype_posteriors[g] * inner_product(component_responsibilities[s][t][k], log_likelihoods[s][g][t][0]);
                    }
                }
            }
        }
        maths::normalise_exp(result[s]);
    }
}

VariationalBayesMixtureMixtureModel::ComponentResponsibilityMatrix
VariationalBayesMixtureMixtureModel::init_responsibilities(const GroupConcentrationVector& group_concentrations,
                                                           const MixtureConcentrationArray& mixture_concentrations,
                                                           const ProbabilityVector& genotype_priors,
                                                           const GroupResponsibilityVector& group_responsibilities,
                                                           const HaplotypeLikelihoodMatrix& log_likelihoods) const
{
    const auto S = log_likelihoods.size();
    const auto T = group_concentrations.size();
    std::size_t K {0};
    for (std::size_t s {0}; s < S; ++s) {
        for (std::size_t t {0}; t < T; ++t) {
            K = std::max(K, mixture_concentrations[s][t].size());
        }
    }
    ComponentResponsibilityMatrix result(S, ComponentResponsibilityVectorArray(T, ComponentResponsibilityVector(K)));
    for (std::size_t s {0}; s < S; ++s) {
        for (std::size_t t {0}; t < T; ++t) {
            const auto N = log_likelihoods[s][0][0][0].size();
            for (std::size_t k {0}; k < K; ++k) {
                result[s][t][k].resize(N);
            }
        }
    }
    update_responsibilities(result, group_concentrations, mixture_concentrations,
                            genotype_priors, group_responsibilities, log_likelihoods);
    return result;
}

namespace {

auto inner_product(const VariationalBayesMixtureMixtureModel::ProbabilityVector& genotype_posteriors,
                   const VariationalBayesMixtureMixtureModel::GenotypeCombinationLikelihoodVector& log_likelihoods,
                   const std::size_t t, const std::size_t k, const std::size_t n) noexcept
{
    ProbabilityVector::value_type result {0};
    const auto G = genotype_posteriors.size();
    for (std::size_t g {0}; g < G; ++g) {
        result += genotype_posteriors[g] * log_likelihoods[g][t][k][n];
    }
    return result;
}

} // namespace

void
VariationalBayesMixtureMixtureModel::update_responsibilities(ComponentResponsibilityMatrix& result,
                                                             const GroupConcentrationVector& group_concentrations,
                                                             const MixtureConcentrationArray& mixture_concentrations,
                                                             const ProbabilityVector& genotype_posteriors,
                                                             const GroupResponsibilityVector& group_responsibilities,
                                                             const HaplotypeLikelihoodMatrix& log_likelihoods) const
{
    const auto T = group_concentrations.size();
    const auto S = log_likelihoods.size();
    const auto max_K = result[0][0].size();
    for (std::size_t s {0}; s < S; ++s) {
        const auto N = log_likelihoods[s][0][0][0].size();
        for (std::size_t t {0}; t < T; ++t) {
            const auto ln_exp_pi = dirichlet_expectation_log(mixture_concentrations[s][t]);
            const auto K = ln_exp_pi.size();
            for (std::size_t k {0}; k < max_K; ++k) {
                for (std::size_t n {0}; n < N; ++n) {
                    if (t == 0) result[s][t][k][n] = 0;
                    auto w = k < K ? ln_exp_pi[k] + inner_product(genotype_posteriors, log_likelihoods[s], t, k, n) : options_.null_log_probability;
                    result[s][t][k][n] += group_responsibilities[s][t] * w;
                }
            }
            std::vector<double> ln_rho(max_K);
            for (std::size_t n {0}; n < N; ++n) {
                for (std::size_t k {0}; k < max_K; ++k) {
                    ln_rho[k] = result[s][t][k][n];
                }
                const auto ln_rho_norm = log_sum_exp(ln_rho);
                for (std::size_t k {0}; k < max_K; ++k) {
                    result[s][t][k][n] = std::exp(ln_rho[k] - ln_rho_norm);
                }
            }
        }
    }
}

void
VariationalBayesMixtureMixtureModel::update_genotype_log_posteriors(LogProbabilityVector& result,
                                                                    const LogProbabilityVector& genotype_log_priors,
                                                                    const GroupResponsibilityVector& group_responsibilities,
                                                                    const ComponentResponsibilityMatrix& component_responsibilities,
                                                                    const HaplotypeLikelihoodMatrix& log_likelihoods) const
{
    for (std::size_t g {0}; g < genotype_log_priors.size(); ++g) {
        result[g] = genotype_log_priors[g] + marginalise(group_responsibilities, component_responsibilities, log_likelihoods, g);
    }
    maths::normalise_logs(result);
}

VariationalBayesMixtureMixtureModel::LogProbability
VariationalBayesMixtureMixtureModel::marginalise(const GroupResponsibilityVector& group_responsibilities,
                                                 const ComponentResponsibilityMatrix& component_responsibilities,
                                                 const HaplotypeLikelihoodMatrix& log_likelihoods,
                                                 const std::size_t g) const noexcept
{
    const auto T = group_responsibilities.front().size();
    const auto S = component_responsibilities.size();
    LogProbability result {0};
    for (std::size_t s {0}; s < S; ++s) {
        for (std::size_t t {0}; t < T; ++t) {
            result += group_responsibilities[s][t] * marginalise(component_responsibilities[s][t], log_likelihoods[s][g][t]);
        }
    }
    return result;
}

VariationalBayesMixtureMixtureModel::LogProbability
VariationalBayesMixtureMixtureModel::marginalise(const ComponentResponsibilityVector& responsibilities,
                                                 const HaplotypeLikelihoodVector& log_likelihoods) const noexcept
{
    const auto K = log_likelihoods.size();
    const auto max_K = responsibilities.size();
    assert(K <= max_K);
    LogProbability result {0};
    for (std::size_t k {0}; k < max_K; ++k) {
        if (k < K) {
            result += inner_product(responsibilities[k], log_likelihoods[k]);
        } else {
            result += inner_product(responsibilities[k], log_likelihoods[0]);
        }
    }
    return result;
}

void
VariationalBayesMixtureMixtureModel::update_group_concentrations(GroupConcentrationVector& result,
                                                                 const GroupConcentrationVector& prior_group_concentrations,
                                                                 const GroupResponsibilityVector& group_responsibilities) const
{
    assert(result.size() == prior_group_concentrations.size());
    const auto T = prior_group_concentrations.size();
    const auto S = group_responsibilities.size();
    for (std::size_t t {0}; t < T; ++t) {
        result[t] = prior_group_concentrations[t];
        for (std::size_t s {0}; s < S; ++s) {
            result[t] += group_responsibilities[s][t];
        }
    }
}

void
VariationalBayesMixtureMixtureModel::update_mixture_concentrations(MixtureConcentrationArray& result,
                                                                   const MixtureConcentrationArray& prior_mixture_concentrations,
                                                                   const GroupResponsibilityVector& group_responsibilities,
                                                                   const ComponentResponsibilityMatrix& component_responsibilities) const
{
    const auto S = prior_mixture_concentrations.size();
    const auto T = group_responsibilities.front().size();
    for (std::size_t s {0}; s < S; ++s) {
        for (std::size_t t {0}; t < T; ++t) {
            const auto K = prior_mixture_concentrations[s][t].size();
            for (std::size_t k {0}; k < K; ++k) {
                result[s][t][k] = prior_mixture_concentrations[s][t][k] + group_responsibilities[s][t] * sum(component_responsibilities[s][t][k]);
            }
        }
    }
}

namespace {

template <typename T>
auto shannon_entropy(const std::vector<T>& probabilities) noexcept
{
    return -std::accumulate(std::cbegin(probabilities), std::cend(probabilities), T {0},
                            [] (const auto curr, const auto p) noexcept { return curr + (p > 0 ? p * std::log(p) : 0.0); });
}

template <typename T>
auto shannon_entropy(const std::vector<std::vector<T>>& probabilities) noexcept
{
    return std::accumulate(std::cbegin(probabilities), std::cend(probabilities), T {0},
                           [] (const auto curr, const auto& ps) noexcept { return curr + shannon_entropy(ps); });
}

} // namespace

double
VariationalBayesMixtureMixtureModel::calculate_evidence(const GroupConcentrationVector& prior_group_concentrations,
                                                        const GroupConcentrationVector& posterior_group_concentrations,
                                                        const MixtureConcentrationArray& prior_mixture_concentrations,
                                                        const MixtureConcentrationArray& posterior_mixture_concentrations,
                                                        const LogProbabilityVector& genotype_log_priors,
                                                        const LogProbabilityVector& genotype_log_posteriors,
                                                        const ProbabilityVector& genotype_posteriors,
                                                        const GroupResponsibilityVector& group_responsibilities,
                                                        const ComponentResponsibilityMatrix& component_responsibilities,
                                                        const HaplotypeLikelihoodMatrix& log_likelihoods) const
{
    const auto G = genotype_log_priors.size();
    const auto S = prior_mixture_concentrations.size();
    const auto T = group_responsibilities.front().size();
    double result {0};
    for (std::size_t g {0}; g < G; ++g) {
        auto w = genotype_log_priors[g] - genotype_log_posteriors[g];
        for (std::size_t s {0}; s < S; ++s) {
            double ss {0};
            for (std::size_t t {0}; t < T; ++t) {
                ss += group_responsibilities[s][t] * marginalise(component_responsibilities[s][t], log_likelihoods[s][g][t]);
            }
            w += ss;
        }
        result += genotype_posteriors[g] * w;
    }    
    for (std::size_t s {0}; s < S; ++s) {
        result += shannon_entropy(group_responsibilities[s]);
        for (std::size_t t {0}; t < T; ++t) {
            result += group_responsibilities[s][t] * shannon_entropy(component_responsibilities[s][t]);
            result += maths::log_beta(posterior_mixture_concentrations[s][t]) - maths::log_beta(prior_mixture_concentrations[s][t]);
        }
    }    
    result += maths::log_beta(posterior_group_concentrations) - maths::log_beta(prior_group_concentrations);    
    return result;
}

namespace {

template <typename Range>
void check_normalisation(Range& probabilities) noexcept
{
    const auto mass = std::accumulate(std::cbegin(probabilities), std::cend(probabilities), 0.0);
    if (mass > 1.0) for (auto& p : probabilities) p /= mass;
}

} // namespace

void
VariationalBayesMixtureMixtureModel::compute_evidence_weighted_latents(ProbabilityVector& genotype_posteriors,
                                                                       GroupResponsibilityVector& group_responsibilities,
                                                                       const std::vector<PointInferences>& latents) const
{
    assert(!latents.empty());
    const auto num_genotypes = latents.front().latents.genotype_posteriors.size();
    std::vector<std::size_t> map_genotypes(num_genotypes, latents.size());
    for (std::size_t i {0}; i < latents.size(); ++i) {
        const auto map_genotype_idx = max_element_index(latents[i].latents.genotype_posteriors);
        if (map_genotypes[map_genotype_idx] == latents.size()
         || latents[i].approx_log_evidence > latents[map_genotypes[map_genotype_idx]].approx_log_evidence) {
            map_genotypes[map_genotype_idx] = i;
        }
    }
    std::vector<std::size_t> modes {};
    modes.reserve(latents.size());
    for (std::size_t g {0}; g < num_genotypes; ++g) {
        if (map_genotypes[g] < latents.size()) {
            modes.push_back(map_genotypes[g]);
        }
    }
    std::vector<double> mode_weights(modes.size());
    std::transform(std::cbegin(modes), std::cend(modes), std::begin(mode_weights),
                   [&] (auto mode) { return latents[mode].approx_log_evidence; });
    maths::normalise_exp(mode_weights);
    genotype_posteriors.resize(num_genotypes);
    const auto num_samples = latents.front().latents.group_responsibilities.size();
    const auto num_groups = latents.front().latents.group_responsibilities.front().size();
    group_responsibilities.resize(num_samples, Sigma(num_groups));
    for (std::size_t i {0}; i < modes.size(); ++i) {
        const auto mode_weight = mode_weights[i];
        const auto& mode_genotype_posteriors = latents[modes[i]].latents.genotype_posteriors;
        std::transform(std::cbegin(mode_genotype_posteriors), std::cend(mode_genotype_posteriors),
                       std::cbegin(genotype_posteriors), std::begin(genotype_posteriors),
                       [mode_weight] (auto seed_posterior, auto curr_posterior) {
                           return curr_posterior + mode_weight * seed_posterior;
                       });
        const auto& mode_group_responsabilities = latents[modes[i]].latents.group_responsibilities;
        for (std::size_t s {0}; s < num_samples; ++s) {
            for (std::size_t t {0}; t < num_groups; ++t) {
                group_responsibilities[s][t] += mode_weight * mode_group_responsabilities[s][t];
            }
        }
    }
    check_normalisation(genotype_posteriors);
    for (auto& responsibilities : group_responsibilities)  check_normalisation(responsibilities);
}

void VariationalBayesMixtureMixtureModel::print_concentrations(const GroupConcentrationVector& concentrations) const
{
    for (std::size_t t {0}; t < concentrations.size(); ++t) {
        std::cout << "t: " << t << " = " << concentrations[t] << std::endl;
    }
}

void VariationalBayesMixtureMixtureModel::print_concentrations(const MixtureConcentrationArray& concentrations) const
{
    for (std::size_t s {0}; s < concentrations.size(); ++s) {
        for (std::size_t t {0}; t < concentrations[s].size(); ++t) {
            for (std::size_t k {0}; k < concentrations[s][t].size(); ++k) {
                std::cout << "s: " << s << " t: " << t << " k: " << k << " = " << concentrations[s][t][k] << std::endl;
            }
        }
    }
}

void VariationalBayesMixtureMixtureModel::print(const ProbabilityVector& probabilities) const
{
    for (std::size_t g {0}; g < probabilities.size(); ++g) {
        std::cout << "g: " << g << " = " << probabilities[g] << std::endl;
    }
}

void VariationalBayesMixtureMixtureModel::print(const GroupResponsibilityVector& responsibilities) const
{
    for (std::size_t s {0}; s < responsibilities.size(); ++s) {
        for (std::size_t t {0}; t < responsibilities[s].size(); ++t) {
            std::cout << "s: " << s << " t: " << t << " = " << responsibilities[s][t] << std::endl;
        }
    }
}

void VariationalBayesMixtureMixtureModel::print(const ComponentResponsibilityMatrix& responsibilities) const
{
    for (std::size_t s {0}; s < responsibilities.size(); ++s) {
        for (std::size_t t {0}; t < responsibilities[s].size(); ++t) {
            for (std::size_t k {0}; k < responsibilities[s].size(); ++k) {
                for (std::size_t n {0}; n < responsibilities[s][k].size(); ++n) {
                    std::cout << "s: " << s << " t: " << t <<  " k: " << k << " n: " << n << " = " << responsibilities[s][t][k][n] << std::endl;
                }
            }
        }
    }
}

void VariationalBayesMixtureMixtureModel::print(const HaplotypeLikelihoodMatrix& log_likelihoods) const
{
    for (std::size_t s {0}; s < log_likelihoods.size(); ++s) {
        for (std::size_t g {0}; g < log_likelihoods[s].size(); ++g) {
            for (std::size_t t {0}; t < log_likelihoods[s][g].size(); ++t) {
                for (std::size_t k {0}; k < log_likelihoods[s][g][t].size(); ++k) {
                    for (std::size_t n {0}; n < log_likelihoods[s][g][t][k].size(); ++n) {
                        std::cout << "s: " << s <<  " g: " << g << " t: " << t << " k: " << k
                                  << " n: " << n << " = " << log_likelihoods[s][g][t][k][n] << std::endl;
                    }
                }
            }
        }
    }
}

} // namespace model
} // namespace octopus
