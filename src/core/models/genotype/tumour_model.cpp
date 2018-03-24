// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "tumour_model.hpp"

#include <algorithm>
#include <numeric>
#include <iterator>
#include <functional>
#include <cstddef>
#include <cmath>
#include <cassert>

#include "utils/maths.hpp"
#include "logging/logging.hpp"
#include "germline_likelihood_model.hpp"
#include "variational_bayes_mixture_model.hpp"

namespace octopus { namespace model {

TumourModel::TumourModel(std::vector<SampleName> samples, Priors priors)
: TumourModel {std::move(samples), std::move(priors), AlgorithmParameters {}}
{}

TumourModel::TumourModel(std::vector<SampleName> samples, Priors priors, AlgorithmParameters parameters)
: samples_ {std::move(samples)}
, priors_ {std::move(priors)}
, parameters_ {parameters}
{}

const TumourModel::Priors& TumourModel::priors() const noexcept
{
    return priors_;
}

namespace {

template <std::size_t K>
TumourModel::InferredLatents
run_variational_bayes(const std::vector<SampleName>& samples,
                      const std::vector<CancerGenotype<Haplotype>>& genotypes,
                      const TumourModel::Priors& priors,
                      const HaplotypeLikelihoodCache& haplotype_log_likelihoods,
                      const VariationalBayesParameters& params);

template <std::size_t K>
TumourModel::InferredLatents
run_variational_bayes(const std::vector<SampleName>& samples,
                      const std::vector<CancerGenotype<Haplotype>>& genotypes,
                      const std::vector<std::pair<std::vector<unsigned>, unsigned>>& genotype_indices,
                      const TumourModel::Priors& priors,
                      const HaplotypeLikelihoodCache& haplotype_log_likelihoods,
                      const VariationalBayesParameters& params);

} // namespace

TumourModel::InferredLatents
TumourModel::evaluate(const std::vector<CancerGenotype<Haplotype>>& genotypes,
                      const HaplotypeLikelihoodCache& haplotype_likelihoods) const
{
    assert(!genotypes.empty());
    const VariationalBayesParameters vb_params {parameters_.epsilon, parameters_.max_iterations};
    auto ploidy = genotypes.front().ploidy();
    assert(ploidy < 3);
    if (ploidy == 1) {
        return run_variational_bayes<2>(samples_, genotypes, priors_, haplotype_likelihoods, vb_params);
    }
    return run_variational_bayes<3>(samples_, genotypes, priors_, haplotype_likelihoods, vb_params);
}

TumourModel::InferredLatents
TumourModel::evaluate(const std::vector<CancerGenotype<Haplotype>>& genotypes,
                      const std::vector<std::pair<std::vector<unsigned>, unsigned>>& genotype_indices,
                      const HaplotypeLikelihoodCache& haplotype_likelihoods) const
{
    assert(!genotypes.empty());
    assert(genotypes.size() == genotype_indices.size());
    const VariationalBayesParameters vb_params {parameters_.epsilon, parameters_.max_iterations};
    auto ploidy = genotypes.front().ploidy();
    assert(ploidy < 3);
    if (ploidy == 1) {
        return run_variational_bayes<2>(samples_, genotypes, genotype_indices, priors_, haplotype_likelihoods, vb_params);
    }
    return run_variational_bayes<3>(samples_, genotypes, genotype_indices, priors_, haplotype_likelihoods, vb_params);
}

namespace {

template <std::size_t K>
VBAlpha<K> flatten(const TumourModel::Priors::GenotypeMixturesDirichletAlphas& alpha)
{
    VBAlpha<K> result {};
    std::copy_n(std::cbegin(alpha), K, std::begin(result));
    return result;
}

template <std::size_t K>
VBAlphaVector<K> flatten(const TumourModel::Priors::GenotypeMixturesDirichletAlphaMap& alphas,
                         const std::vector<SampleName>& samples)
{
    VBAlphaVector<K> result(samples.size());
    std::transform(std::cbegin(samples), std::cend(samples), std::begin(result),
                   [&alphas] (const auto& sample) { return flatten<K>(alphas.at(sample)); });
    return result;
}

template <std::size_t K>
VBGenotype<K>
flatten(const CancerGenotype<Haplotype>& genotype, const SampleName& sample,
        const HaplotypeLikelihoodCache& haplotype_likelihoods)
{
    VBGenotype<K> result {};
    const Genotype<Haplotype>& germline_genotype{genotype.germline_genotype()};
    assert(germline_genotype.ploidy() == (K - 1));
    std::transform(std::cbegin(germline_genotype), std::cend(germline_genotype), std::begin(result),
                   [&sample, &haplotype_likelihoods] (const Haplotype& haplotype)
                   -> std::reference_wrapper<const VBReadLikelihoodArray::BaseType> {
                       return std::cref(haplotype_likelihoods(sample, haplotype));
                   });
    result.back() = haplotype_likelihoods(sample, genotype.somatic_element());
    return result;
}

template <std::size_t K>
VBGenotypeVector<K>
flatten(const std::vector<CancerGenotype<Haplotype>>& genotypes, const SampleName& sample,
        const HaplotypeLikelihoodCache& haplotype_likelihoods)
{
    VBGenotypeVector<K> result(genotypes.size());
    std::transform(std::cbegin(genotypes), std::cend(genotypes), std::begin(result),
                   [&sample, &haplotype_likelihoods] (const auto& genotype) {
                       return flatten<K>(genotype, sample, haplotype_likelihoods);
                   });
    return result;
}

template <std::size_t K>
VBReadLikelihoodMatrix<K>
flatten(const std::vector<CancerGenotype<Haplotype>>& genotypes,
        const std::vector<SampleName>& samples,
        const HaplotypeLikelihoodCache& haplotype_likelihoods)
{
    VBReadLikelihoodMatrix<K> result {};
    result.reserve(samples.size());
    std::transform(std::cbegin(samples), std::cend(samples), std::back_inserter(result),
                   [&genotypes, &haplotype_likelihoods] (const auto& sample) {
                       return flatten<K>(genotypes, sample, haplotype_likelihoods);
                   });
    return result;
}
template <std::size_t K>
TumourModel::Latents::GenotypeMixturesDirichletAlphas expand(VBAlpha<K>& alpha)
{
    return TumourModel::Latents::GenotypeMixturesDirichletAlphas(std::begin(alpha), std::end(alpha));
}

template <std::size_t K>
TumourModel::Latents::GenotypeMixturesDirichletAlphaMap
expand(const std::vector<SampleName>& samples, VBAlphaVector<K>&& alphas)
{
    TumourModel::Latents::GenotypeMixturesDirichletAlphaMap result {};
    std::transform(std::cbegin(samples), std::cend(samples), std::begin(alphas),
                   std::inserter(result, std::begin(result)),
                   [] (const auto& sample, auto&& vb_alpha) {
                       return std::make_pair(sample, expand(vb_alpha));
                   });
    return result;
}

template <std::size_t K>
TumourModel::InferredLatents
expand(const std::vector<SampleName>& samples, VBLatents<K>&& inferred_latents, double evidence)
{
    TumourModel::Latents posterior_latents {
        std::move(inferred_latents.genotype_posteriors),
        expand(samples, std::move(inferred_latents.alphas))
    };
    return {std::move(posterior_latents), evidence};
}

auto calculate_log_priors(const std::vector<std::pair<std::vector<unsigned>, unsigned>>& genotype_indices,
                          const CancerGenotypePriorModel& model)
{
    
    std::vector<double> result(genotype_indices.size());
    std::transform(std::cbegin(genotype_indices), std::cend(genotype_indices), std::begin(result),
                   [&] (const auto& p) { return model.evaluate(p.first, p.second); });
    maths::normalise_logs(result);
    return result;
}

auto compute_germline_log_likelihoods(const SampleName& sample,
                                      const std::vector<CancerGenotype<Haplotype>>& genotypes,
                                      const HaplotypeLikelihoodCache& haplotype_log_likelihoods)
{
    haplotype_log_likelihoods.prime(sample);
    const GermlineLikelihoodModel likelihood_model {haplotype_log_likelihoods};
    std::vector<double> result(genotypes.size());
    std::transform(std::cbegin(genotypes), std::cend(genotypes), std::begin(result),
                   [&] (const auto& genotype) { return likelihood_model.evaluate(genotype.germline_genotype()); });
    return result;
}

auto compute_demoted_log_likelihoods(const SampleName& sample,
                                     const std::vector<CancerGenotype<Haplotype>>& genotypes,
                                     const HaplotypeLikelihoodCache& haplotype_log_likelihoods)
{
    assert(!genotypes.empty());
    haplotype_log_likelihoods.prime(sample);
    const GermlineLikelihoodModel likelihood_model {haplotype_log_likelihoods};
    std::vector<double> result(genotypes.size());
    std::transform(std::cbegin(genotypes), std::cend(genotypes), std::begin(result),
                   [&] (const auto& genotype) { return likelihood_model.evaluate(demote(genotype)); });
    return result;
}

auto compute_log_posteriors(const LogProbabilityVector& log_priors, const LogProbabilityVector& log_likelihoods)
{
    assert(log_priors.size() == log_likelihoods.size());
    LogProbabilityVector result(log_priors.size());
    std::transform(std::cbegin(log_priors), std::cend(log_priors), std::cbegin(log_likelihoods), std::begin(result),
                   [] (auto prior, auto likelihood) noexcept { return prior + likelihood; });
    maths::normalise_logs(result);
    return result;
}

LogProbabilityVector log_uniform_dist(const std::size_t n)
{
    return LogProbabilityVector(n, -std::log(static_cast<double>(n)));
}

auto make_point_seed(const std::size_t num_genotypes, const std::size_t n, const double p = 0.99)
{
    LogProbabilityVector result(num_genotypes, std::log((1 - p) / (num_genotypes - 1)));
    result[n] = std::log(p);
    return result;
}

auto make_range_seed(const std::size_t num_genotypes, const std::size_t begin, const std::size_t n, const double p = 0.99)
{
    LogProbabilityVector result(num_genotypes, std::log((1 - p) / (num_genotypes - n)));
    std::fill_n(std::next(std::begin(result), begin), n, std::log(p / n));
    return result;
}

auto generate_seeds(const std::vector<SampleName>& samples,
                    const std::vector<CancerGenotype<Haplotype>>& genotypes,
                    const LogProbabilityVector& genotype_log_priors,
                    const HaplotypeLikelihoodCache& haplotype_log_likelihoods)
{
    std::vector<LogProbabilityVector> result {};
    result.reserve(2 + 3 * samples.size());
    result.push_back(genotype_log_priors);
    result.push_back(log_uniform_dist(genotypes.size()));
    for (const auto& sample : samples) {
        auto log_likelihoods = compute_germline_log_likelihoods(sample, genotypes, haplotype_log_likelihoods);
        result.push_back(compute_log_posteriors(genotype_log_priors, log_likelihoods));
        maths::normalise_logs(log_likelihoods);
        result.push_back(std::move(log_likelihoods));
        auto demoted_log_likelihoods = compute_demoted_log_likelihoods(sample, genotypes, haplotype_log_likelihoods);
        result.push_back(compute_log_posteriors(genotype_log_priors, demoted_log_likelihoods));
        maths::normalise_logs(demoted_log_likelihoods);
        result.push_back(std::move(demoted_log_likelihoods));
    }
    return result;
}

template <std::size_t K>
TumourModel::InferredLatents
run_variational_bayes(const std::vector<SampleName>& samples,
                      const std::vector<CancerGenotype<Haplotype>>& genotypes,
                      const TumourModel::Priors::GenotypeMixturesDirichletAlphaMap& prior_alphas,
                      const std::vector<double>& genotype_log_priors,
                      const HaplotypeLikelihoodCache& haplotype_log_likelihoods,
                      const VariationalBayesParameters& params,
                      std::vector<std::vector<double>>&& seeds)
{
    const auto vb_prior_alphas = flatten<K>(prior_alphas, samples);
    const auto log_likelihoods = flatten<K>(genotypes, samples, haplotype_log_likelihoods);
    auto p = run_variational_bayes(vb_prior_alphas, genotype_log_priors, log_likelihoods, params, std::move(seeds));
    return expand(samples, std::move(p.first), p.second);
}

// Main entry point

template <std::size_t K>
TumourModel::InferredLatents
run_variational_bayes(const std::vector<SampleName>& samples,
                      const std::vector<CancerGenotype<Haplotype>>& genotypes,
                      const TumourModel::Priors& priors,
                      const HaplotypeLikelihoodCache& haplotype_log_likelihoods,
                      const VariationalBayesParameters& params)
{
    const auto genotype_log_priors = calculate_log_priors(genotypes, priors.genotype_prior_model);
    auto seeds = generate_seeds(samples, genotypes, genotype_log_priors, haplotype_log_likelihoods);
    return run_variational_bayes<K>(samples, genotypes, priors.alphas, genotype_log_priors,
                                    haplotype_log_likelihoods, params, std::move(seeds));
}

template <std::size_t K>
TumourModel::InferredLatents
run_variational_bayes(const std::vector<SampleName>& samples,
                      const std::vector<CancerGenotype<Haplotype>>& genotypes,
                      const std::vector<std::pair<std::vector<unsigned>, unsigned>>& genotype_indices,
                      const TumourModel::Priors& priors,
                      const HaplotypeLikelihoodCache& haplotype_log_likelihoods,
                      const VariationalBayesParameters& params)
{
    const auto genotype_log_priors = calculate_log_priors(genotype_indices, priors.genotype_prior_model);
    auto seeds = generate_seeds(samples, genotypes, genotype_log_priors, haplotype_log_likelihoods);
    return run_variational_bayes<K>(samples, genotypes, priors.alphas, genotype_log_priors,
                                    haplotype_log_likelihoods, params, std::move(seeds));
}

} // namespace

} // namespace model
} // namespace octopus
