// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "variational_bayes_mixture_mixture_model.hpp"

#include <utility>
#include <cmath>
#include <iterator>
#include <algorithm>
#include <numeric>
#include <iostream>

#include <boost/math/special_functions/digamma.hpp>

#include "utils/maths.hpp"

namespace octopus { namespace model {

VariationalBayesMixtureMixtureModel::VariationalBayesMixtureMixtureModel(Options options)
: options_ {std::move(options)}
{}

VariationalBayesMixtureMixtureModel::Inferences
VariationalBayesMixtureMixtureModel::evaluate(const LogProbabilityVector& genotype_log_priors,
                                              const HaplotypeLikelihoodMatrix& log_likelihoods,
                                              const GroupConcentrationVector& group_concentrations,
                                              const MixtureConcentrationArray& mixture_concentrations,
                                              std::vector<LogProbabilityVector> seeds) const
{
    Inferences result {};
    result.approx_log_evidence = std::numeric_limits<double>::lowest();
    for (auto& seed : seeds) {
        auto inferences = evaluate(genotype_log_priors, log_likelihoods, group_concentrations, mixture_concentrations, seed);
        if (inferences.approx_log_evidence > result.approx_log_evidence) {
            result = std::move(inferences);
        }
    }
    return result;
}

VariationalBayesMixtureMixtureModel::Inferences
VariationalBayesMixtureMixtureModel::evaluate(const LogProbabilityVector& genotype_log_priors,
                                              const HaplotypeLikelihoodMatrix& log_likelihoods,
                                              const double group_concentration, const double mixture_concentration,
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
    return evaluate(genotype_log_priors, log_likelihoods, group_concentrations, mixture_concentrations, std::move(seeds));
}

// Private methods

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
T log_sum_exp(const std::vector<T>& logs)
{
    return maths::log_sum_exp(logs);
}

} // namespace

VariationalBayesMixtureMixtureModel::Inferences
VariationalBayesMixtureMixtureModel::evaluate(const LogProbabilityVector& genotype_log_priors,
                                              const HaplotypeLikelihoodMatrix& log_likelihoods,
                                              const GroupConcentrationVector& prior_group_concentrations,
                                              const MixtureConcentrationArray& prior_mixture_concentrations,
                                              LogProbabilityVector& genotype_log_posteriors) const
{
    auto genotype_posteriors = exp(genotype_log_posteriors);
    auto posterior_group_concentrations = prior_group_concentrations;
    auto posterior_mixture_concentrations = prior_mixture_concentrations;
    auto group_responsabilities = init_responsabilities(prior_group_concentrations, prior_mixture_concentrations,
                                                        genotype_posteriors, log_likelihoods);
    auto component_responsabilities = init_responsabilities(prior_group_concentrations, prior_mixture_concentrations,
                                                            genotype_posteriors, group_responsabilities, log_likelihoods);
    auto prev_evidence = std::numeric_limits<double>::lowest();
    for (unsigned i {0}; i < options_.max_iterations; ++i) {
        update_genotype_log_posteriors(genotype_log_posteriors, genotype_log_priors,
                                       group_responsabilities, component_responsabilities,
                                       log_likelihoods);
        exp(genotype_log_posteriors, genotype_posteriors);
        update_group_concentrations(posterior_group_concentrations, prior_group_concentrations, group_responsabilities);
        update_mixture_concentrations(posterior_mixture_concentrations, prior_mixture_concentrations,
                                      group_responsabilities, component_responsabilities);
        auto curr_evidence = calculate_evidence(prior_group_concentrations, posterior_group_concentrations,
                                                prior_mixture_concentrations, posterior_mixture_concentrations,
                                                genotype_log_priors, genotype_log_posteriors, genotype_posteriors,
                                                group_responsabilities, component_responsabilities,
                                                log_likelihoods);
        //assert(curr_evidence + 2 * options_.epsilon >= prev_evidence);
        if (curr_evidence <= prev_evidence || (curr_evidence - prev_evidence) < options_.epsilon) {
            prev_evidence = curr_evidence;
            break;
        }
        prev_evidence = curr_evidence;
        update_responsabilities(group_responsabilities, posterior_group_concentrations, posterior_mixture_concentrations,
                                genotype_posteriors, component_responsabilities, log_likelihoods);
        update_responsabilities(component_responsabilities, posterior_group_concentrations, posterior_mixture_concentrations,
                                genotype_posteriors, group_responsabilities, log_likelihoods);
    }
    Inferences result {};
    result.genotype_log_posteriors = std::move(genotype_log_posteriors);
    result.genotype_posteriors = std::move(genotype_posteriors);
    result.group_concentrations = std::move(posterior_group_concentrations);
    result.mixture_concentrations = std::move(posterior_mixture_concentrations);
    result.group_responsabilities = std::move(group_responsabilities);
    result.approx_log_evidence = prev_evidence;
    return result;
}

VariationalBayesMixtureMixtureModel::GroupResponsabilityVector
VariationalBayesMixtureMixtureModel::init_responsabilities(const GroupConcentrationVector& group_concentrations,
                                                           const MixtureConcentrationArray& mixture_concentrations,
                                                           const ProbabilityVector& genotype_priors,
                                                           const HaplotypeLikelihoodMatrix& log_likelihoods) const
{
    const auto S = log_likelihoods.size();
    const auto T = group_concentrations.size();
    GroupResponsabilityVector result(S, Sigma(T));
    const auto G = genotype_priors.size();
    const auto beta_0 = sum(group_concentrations);
    std::vector<GroupConcentrationVector::value_type> ln_ex_psi(T);
    std::transform(std::cbegin(group_concentrations), std::cend(group_concentrations), std::begin(ln_ex_psi),
                   [=] (auto beta) { return digamma_diff(beta, beta_0); });
    std::size_t max_K {0};
    for (std::size_t s {0}; s < S; ++s) {
        for (std::size_t t {0}; t < T; ++t) {
            max_K = std::max(max_K, mixture_concentrations[s][t].size());
        }
    }
    for (std::size_t s {0}; s < S; ++s) {
        const auto N = log_likelihoods[s][0][0][0].size();
        const auto tau = 1.0 / max_K;
        const auto tau_sum = N * tau;
        for (std::size_t t {0}; t < T; ++t) {
            result[s][t] = ln_ex_psi[t];
            const auto K = mixture_concentrations[s][t].size();
            std::vector<double> ln_ex_pi(K);
            const auto alpha_0 = sum(mixture_concentrations[s][t]);
            std::transform(std::cbegin(mixture_concentrations[s][t]), std::cend(mixture_concentrations[s][t]),
                           std::begin(ln_ex_pi), [=] (auto alpha) { return digamma_diff(alpha, alpha_0); });
            for (std::size_t k {0}; k < K; ++k) {
                result[s][t] += ln_ex_pi[k] * tau_sum;
                for (std::size_t g {0}; g < G; ++g) {
                    result[s][t] += genotype_priors[g] * tau * sum(log_likelihoods[s][g][t][k]);
                }
            }
        }
        maths::normalise_exp(result[s]);
    }
    return result;
}

void
VariationalBayesMixtureMixtureModel::update_responsabilities(GroupResponsabilityVector& result,
                                                             const GroupConcentrationVector& group_concentrations,
                                                             const MixtureConcentrationArray& mixture_concentrations,
                                                             const ProbabilityVector& genotype_posteriors,
                                                             const ComponentResponsabilityMatrix& component_responsabilities,
                                                             const HaplotypeLikelihoodMatrix& log_likelihoods) const
{
    const auto T = group_concentrations.size();
    const auto S = log_likelihoods.size();
    const auto G = genotype_posteriors.size();
    const auto beta_0 = sum(group_concentrations);
    std::vector<GroupConcentrationVector::value_type> ln_ex_psi(T);
    std::transform(std::cbegin(group_concentrations), std::cend(group_concentrations), std::begin(ln_ex_psi),
                   [=] (auto beta) { return digamma_diff(beta, beta_0); });
    for (std::size_t s {0}; s < S; ++s) {
        std::vector<double> component_responsability_sums(component_responsabilities[s].size());
        std::transform(std::cbegin(component_responsabilities[s]), std::cend(component_responsabilities[s]),
                       std::begin(component_responsability_sums), [] (const auto& taus) { return sum(taus); });
        for (std::size_t t {0}; t < T; ++t) {
            result[s][t] = ln_ex_psi[t];
            const auto K = mixture_concentrations[s][t].size();
            std::vector<double> ln_ex_pi(K);
            const auto alpha_0 = sum(mixture_concentrations[s][t]);
            std::transform(std::cbegin(mixture_concentrations[s][t]), std::cend(mixture_concentrations[s][t]), std::begin(ln_ex_pi),
                           [=] (auto alpha) { return digamma_diff(alpha, alpha_0); });
            for (std::size_t k {0}; k < K; ++k) {
                result[s][t] += ln_ex_pi[k] * component_responsability_sums[k];
                for (std::size_t g {0}; g < G; ++g) {
                    result[s][t] += genotype_posteriors[g] * inner_product(component_responsabilities[s][k], log_likelihoods[s][g][t][k]);
                }
            }
        }
        maths::normalise_exp(result[s]);
    }
}

VariationalBayesMixtureMixtureModel::ComponentResponsabilityMatrix
VariationalBayesMixtureMixtureModel::init_responsabilities(const GroupConcentrationVector& group_concentrations,
                                                           const MixtureConcentrationArray& mixture_concentrations,
                                                           const ProbabilityVector& genotype_priors,
                                                           const GroupResponsabilityVector& group_responsabilities,
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
    ComponentResponsabilityMatrix result(S, ComponentResponsabilityVector(K));
    for (std::size_t s {0}; s < S; ++s) {
        const auto N = log_likelihoods[s][0][0][0].size();
        for (std::size_t k {0}; k < K; ++k) {
            result[s][k].resize(N);
        }
    }
    update_responsabilities(result, group_concentrations, mixture_concentrations,
                            genotype_priors, group_responsabilities, log_likelihoods);
    return result;
}

void
VariationalBayesMixtureMixtureModel::update_responsabilities(ComponentResponsabilityMatrix& result,
                                                             const GroupConcentrationVector& group_concentrations,
                                                             const MixtureConcentrationArray& mixture_concentrations,
                                                             const ProbabilityVector& genotype_posteriors,
                                                             const GroupResponsabilityVector& group_responsabilities,
                                                             const HaplotypeLikelihoodMatrix& log_likelihoods) const
{
    const auto T = group_concentrations.size();
    const auto S = log_likelihoods.size();
    const auto G = genotype_posteriors.size();
    const auto max_K = result[0].size();
    for (std::size_t s {0}; s < S; ++s) {
        const auto N = log_likelihoods[s][0][0][0].size();
        for (std::size_t t {0}; t < T; ++t) {
            const auto K = mixture_concentrations[s][t].size();
            std::vector<double> ln_exp_pi(K);
            const auto alpha_0 = sum(mixture_concentrations[s][t]);
            std::transform(std::cbegin(mixture_concentrations[s][t]), std::cend(mixture_concentrations[s][t]), std::begin(ln_exp_pi),
                           [=] (auto alpha) { return digamma_diff(alpha, alpha_0); });
            for (std::size_t n {0}; n < N; ++n) {
                for (std::size_t k {0}; k < K; ++k) {
                    auto tmp = ln_exp_pi[k];
                    for (std::size_t g {0}; g < G; ++g) {
                        tmp += genotype_posteriors[g] * log_likelihoods[s][g][t][k][n];
                    }
                    result[s][k][n] += group_responsabilities[s][t] * tmp;
                }
            }
        }
        std::vector<double> ln_rho(max_K);
        for (std::size_t n {0}; n < N; ++n) {
            for (std::size_t k {0}; k < max_K; ++k) {
                ln_rho[k] = result[s][k][n];
            }
            const auto ln_rho_norm = log_sum_exp(ln_rho);
            for (std::size_t k {0}; k < max_K; ++k) {
                result[s][k][n] = std::exp(ln_rho[k] - ln_rho_norm);
            }
        }
    }
}

void
VariationalBayesMixtureMixtureModel::update_genotype_log_posteriors(LogProbabilityVector& result,
                                                                    const LogProbabilityVector& genotype_log_priors,
                                                                    const GroupResponsabilityVector& group_responsabilities,
                                                                    const ComponentResponsabilityMatrix& component_responsabilities,
                                                                    const HaplotypeLikelihoodMatrix& log_likelihoods) const
{
    for (std::size_t g {0}; g < genotype_log_priors.size(); ++g) {
        result[g] = genotype_log_priors[g] + marginalise(group_responsabilities, component_responsabilities, log_likelihoods, g);
    }
    maths::normalise_logs(result);
}

VariationalBayesMixtureMixtureModel::LogProbability
VariationalBayesMixtureMixtureModel::marginalise(const GroupResponsabilityVector& group_responsabilities,
                                                 const ComponentResponsabilityMatrix& component_responsabilities,
                                                 const HaplotypeLikelihoodMatrix& log_likelihoods,
                                                 const std::size_t g) const noexcept
{
    const auto T = group_responsabilities.front().size();
    const auto S = component_responsabilities.size();
    LogProbability result {0};
    for (std::size_t s {0}; s < S; ++s) {
        for (std::size_t t {0}; t < T; ++t) {
            result += group_responsabilities[s][t] * marginalise(component_responsabilities[s], log_likelihoods[s][g][t]);
        }
    }
    return result;
}

VariationalBayesMixtureMixtureModel::LogProbability
VariationalBayesMixtureMixtureModel::marginalise(const ComponentResponsabilityVector& responsabilities,
                                                 const HaplotypeLikelihoodVector& log_likelihoods) const noexcept
{
    const auto K = log_likelihoods.size();
    LogProbability result {0};
    for (std::size_t k {0}; k < K; ++k) {
        result += inner_product(responsabilities[k], log_likelihoods[k]);
    }
    return result;
}

void
VariationalBayesMixtureMixtureModel::update_group_concentrations(GroupConcentrationVector& result,
                                                                 const GroupConcentrationVector& prior_group_concentrations,
                                                                 const GroupResponsabilityVector& group_responsabilities) const
{
    assert(result.size() == prior_group_concentrations.size());
    const auto T = prior_group_concentrations.size();
    const auto S = group_responsabilities.size();
    for (std::size_t t {0}; t < T; ++t) {
        result[t] = prior_group_concentrations[t];
        for (std::size_t s {0}; s < S; ++s) {
            result[t] += group_responsabilities[s][t];
        }
    }
}

void
VariationalBayesMixtureMixtureModel::update_mixture_concentrations(MixtureConcentrationArray& result,
                                                                   const MixtureConcentrationArray& prior_mixture_concentrations,
                                                                   const GroupResponsabilityVector& group_responsabilities,
                                                                   const ComponentResponsabilityMatrix& component_responsabilities) const
{
    const auto S = prior_mixture_concentrations.size();
    const auto T = group_responsabilities.front().size();
    for (std::size_t s {0}; s < S; ++s) {
        for (std::size_t t {0}; t < T; ++t) {
            for (std::size_t k {0}; k < prior_mixture_concentrations[s][t].size(); ++k) {
                result[s][t][k] = prior_mixture_concentrations[s][t][k] + group_responsabilities[s][t] * sum(component_responsabilities[s][k]);
            }
        }
    }
}

double
VariationalBayesMixtureMixtureModel::calculate_evidence(const GroupConcentrationVector& prior_group_concentrations,
                                                        const GroupConcentrationVector& posterior_group_concentrations,
                                                        const MixtureConcentrationArray& prior_mixture_concentrations,
                                                        const MixtureConcentrationArray& posterior_mixture_concentrations,
                                                        const LogProbabilityVector& genotype_log_priors,
                                                        const LogProbabilityVector& genotype_log_posteriors,
                                                        const ProbabilityVector& genotype_posteriors,
                                                        const GroupResponsabilityVector& group_responsabilities,
                                                        const ComponentResponsabilityMatrix& component_responsabilities,
                                                        const HaplotypeLikelihoodMatrix& log_likelihoods) const
{
    const auto G = genotype_log_priors.size();
    const auto S = prior_mixture_concentrations.size();
    const auto T = group_responsabilities.front().size();
    const auto max_K = component_responsabilities[0].size();
    double result {0};
    for (std::size_t g {0}; g < G; ++g) {
        auto w = genotype_log_priors[g] - genotype_log_posteriors[g];
        for (std::size_t s {0}; s < S; ++s) {
            const auto N = log_likelihoods[s][0][0][0].size();
            for (std::size_t t {0}; t < T; ++t) {
                auto h = 0;
                const auto K = posterior_mixture_concentrations[s][t].size();
                for (std::size_t k {0}; k < K; ++k) {
                    for (std::size_t n {0}; n < N; ++n) {
                        h += component_responsabilities[s][k][n] * log_likelihoods[s][g][t][k][n];
                    }
                }
                w += group_responsabilities[s][t] * h;
            }
        }
        result += genotype_posteriors[g] * w;
    }
    for (std::size_t s {0}; s < S; ++s) {
        const auto N = log_likelihoods[s][0][0][0].size();
        for (std::size_t n {0}; n < N; ++n) {
            for (std::size_t k {0}; k < max_K; ++k) {
                if (component_responsabilities[s][k][n] > 0) {
                    result -= component_responsabilities[s][k][n] * std::log(component_responsabilities[s][k][n]);
                }
            }
        }
        for (std::size_t t {0}; t < T; ++t) {
            if (group_responsabilities[s][t] > 0) {
                result -= group_responsabilities[s][t] * std::log(group_responsabilities[s][t]);
            }
            result += maths::log_beta(posterior_mixture_concentrations[s][t]) - maths::log_beta(prior_mixture_concentrations[s][t]);
        }
    }
    result += maths::log_beta(posterior_group_concentrations) - maths::log_beta(prior_group_concentrations);
    return result;
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

void VariationalBayesMixtureMixtureModel::print(const GroupResponsabilityVector& responsabilities) const
{
    for (std::size_t s {0}; s < responsabilities.size(); ++s) {
        for (std::size_t t {0}; t < responsabilities[s].size(); ++t) {
            std::cout << "s: " << s << " t: " << t << " = " << responsabilities[s][t] << std::endl;
        }
    }
}

void VariationalBayesMixtureMixtureModel::print(const ComponentResponsabilityMatrix& responsabilities) const
{
    for (std::size_t s {0}; s < responsabilities.size(); ++s) {
        for (std::size_t k {0}; k < responsabilities[s].size(); ++k) {
            for (std::size_t n {0}; n < responsabilities[s][k].size(); ++n) {
                std::cout << "s: " << s <<  " k: " << k << " n: " << n << " = " << responsabilities[s][k][n] << std::endl;
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
