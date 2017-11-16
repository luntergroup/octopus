// Copyright (c) 2017 Daniel Cooke
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

TumourModel::TumourModel(std::vector<SampleName> samples, const unsigned ploidy, Priors priors)
: TumourModel {std::move(samples), ploidy, std::move(priors), AlgorithmParameters {}}
{}

TumourModel::TumourModel(std::vector<SampleName> samples, const unsigned ploidy, Priors priors,
                         AlgorithmParameters parameters)
: samples_ {std::move(samples)}
, ploidy_ {ploidy}
, priors_ {std::move(priors)}
, parameters_ {parameters}
{}

namespace {

template <std::size_t K>
TumourModel::InferredLatents
run_variational_bayes(const std::vector<SampleName>& samples,
                      std::vector<CancerGenotype<Haplotype>>&& genotypes,
                      const TumourModel::Priors& priors,
                      const HaplotypeLikelihoodCache& haplotype_log_likelihoods,
                      const VariationalBayesParameters& params);

} // namespace

TumourModel::InferredLatents
TumourModel::evaluate(std::vector<CancerGenotype<Haplotype>> genotypes,
                      const HaplotypeLikelihoodCache& haplotype_likelihoods) const
{
    assert(!genotypes.empty());
    const VariationalBayesParameters vb_params {parameters_.epsilon, parameters_.max_iterations};
    assert(ploidy_ < 3);
    if (ploidy_ == 1) {
        return run_variational_bayes<2>(samples_, std::move(genotypes), priors_, haplotype_likelihoods, vb_params);
    }
    return run_variational_bayes<3>(samples_, std::move(genotypes), priors_, haplotype_likelihoods, vb_params);
}

namespace {

template <std::size_t K>
CompressedAlpha<K> compress(const TumourModel::Priors::GenotypeMixturesDirichletAlphas& alpha)
{
    CompressedAlpha<K> result;
    std::copy_n(std::cbegin(alpha), K, std::begin(result));
    return result;
}

template <std::size_t K>
CompressedAlphaVector<K> flatten_priors(const TumourModel::Priors& priors, const std::vector<SampleName>& samples)
{
    CompressedAlphaVector<K> result(samples.size());
    std::transform(std::cbegin(samples), std::cend(samples), std::begin(result),
                   [&priors](const auto& sample) {
                       return compress<K>(priors.alphas.at(sample));
                   });
    return result;
}

template <std::size_t K>
CompressedGenotype<K>
compress(const CancerGenotype<Haplotype>& genotype, const SampleName& sample,
         const HaplotypeLikelihoodCache& haplotype_likelihoods)
{
    CompressedGenotype<K> result;
    const Genotype<Haplotype>& germline_genotype{genotype.germline_genotype()};
    assert(germline_genotype.ploidy() == (K - 1));
    std::transform(std::cbegin(germline_genotype), std::cend(germline_genotype),
                   std::begin(result),
                   [&sample, &haplotype_likelihoods](const Haplotype& haplotype)
                   -> std::reference_wrapper<const ReadLikelihoodArray::BaseType> {
                       return std::cref(haplotype_likelihoods(sample, haplotype));
                   });
    result.back() = haplotype_likelihoods(sample, genotype.somatic_element());
    return result;
}

template <std::size_t K>
CompressedGenotypeVector<K>
compress(const std::vector<CancerGenotype<Haplotype>>& genotypes, const SampleName& sample,
         const HaplotypeLikelihoodCache& haplotype_likelihoods)
{
    CompressedGenotypeVector<K> result(genotypes.size());
    std::transform(std::cbegin(genotypes), std::cend(genotypes), std::begin(result),
                   [&sample, &haplotype_likelihoods](const auto& genotype) {
                       return compress<K>(genotype, sample, haplotype_likelihoods);
                   });
    return result;
}

template <std::size_t K>
CompressedReadLikelihoodMatrix<K>
compress(const std::vector<CancerGenotype<Haplotype>>& genotypes,
         const std::vector<SampleName>& samples,
         const HaplotypeLikelihoodCache& haplotype_likelihoods)
{
    CompressedReadLikelihoodMatrix<K> result{};
    result.reserve(samples.size());
    std::transform(std::cbegin(samples), std::cend(samples), std::back_inserter(result),
                   [&genotypes, &haplotype_likelihoods](const auto& sample) {
                       return compress<K>(genotypes, sample, haplotype_likelihoods);
                   });
    return result;
}

TumourModel::Latents::GenotypeProbabilityMap
expand(std::vector<CancerGenotype<Haplotype>>&& genotypes, LogProbabilityVector&& genotype_log_posteriors)
{
    TumourModel::Latents::GenotypeProbabilityMap result{};
    std::transform(std::make_move_iterator(std::begin(genotypes)),
                   std::make_move_iterator(std::end(genotypes)),
                   std::begin(genotype_log_posteriors),
                   std::inserter(result, std::begin(result)),
                   [](auto&& g, auto p) { return std::make_pair(std::move(g), p); });
    return result;
}

template <std::size_t K>
TumourModel::Latents::GenotypeMixturesDirichletAlphas expand(CompressedAlpha<K>& alpha)
{
    return TumourModel::Latents::GenotypeMixturesDirichletAlphas(std::begin(alpha), std::end(alpha));
}

template <std::size_t K>
TumourModel::Latents::GenotypeMixturesDirichletAlphaMap
expand(const std::vector<SampleName>& samples, CompressedAlphaVector<K>&& alphas)
{
    TumourModel::Latents::GenotypeMixturesDirichletAlphaMap result{};
    std::transform(std::cbegin(samples), std::cend(samples), std::begin(alphas),
                   std::inserter(result, std::begin(result)),
                   [](const auto& sample, auto&& compressed_alpha) {
                       return std::make_pair(sample, expand(compressed_alpha));
                   });
    return result;
}

template <std::size_t K>
TumourModel::InferredLatents
expand(const std::vector<SampleName>& samples, std::vector<CancerGenotype<Haplotype>>&& genotypes,
       CompressedLatents<K>&& inferred_latents, const double evidence)
{
    TumourModel::Latents posterior_latents{
    expand(std::move(genotypes), std::move(inferred_latents.genotype_posteriors)),
    expand(samples, std::move(inferred_latents.alphas))
    };
    return TumourModel::InferredLatents {std::move(posterior_latents), evidence};
}

auto compute_germline_log_posteriors(const SampleName& sample,
                                     const std::vector<CancerGenotype<Haplotype>>& genotypes,
                                     const HaplotypeLikelihoodCache& haplotype_log_likelihoods)
{
    assert(!genotypes.empty());
    haplotype_log_likelihoods.prime(sample);
    const GermlineLikelihoodModel likelihood_model{haplotype_log_likelihoods};
    std::vector<double> result(genotypes.size());
    std::transform(std::cbegin(genotypes), std::cend(genotypes), std::begin(result),
                   [&](const auto& genotype) {
                       return likelihood_model.evaluate(genotype.germline_genotype());
                   });
    maths::normalise_logs(result);
    return result;
}

auto compute_germline_log_posteriors(const SampleName& sample,
                                     const std::vector<CancerGenotype<Haplotype>>& genotypes,
                                     const HaplotypeLikelihoodCache& haplotype_log_likelihoods,
                                     const CancerGenotypePriorModel& genotype_prior_model)
{
    assert(!genotypes.empty());
    haplotype_log_likelihoods.prime(sample);
    const GermlineLikelihoodModel likelihood_model{haplotype_log_likelihoods};
    std::vector<double> result(genotypes.size());
    std::transform(std::cbegin(genotypes), std::cend(genotypes), std::begin(result),
                   [&](const auto& genotype) {
                       return genotype_prior_model.evaluate(genotype)
                              + likelihood_model.evaluate(genotype.germline_genotype());
                   });
    maths::normalise_logs(result);
    return result;
}

auto compute_log_posteriors_with_germline_model(const SampleName& sample,
                                                const std::vector<CancerGenotype<Haplotype>>& genotypes,
                                                const HaplotypeLikelihoodCache& haplotype_log_likelihoods)
{
    assert(!genotypes.empty());
    haplotype_log_likelihoods.prime(sample);
    const GermlineLikelihoodModel likelihood_model{haplotype_log_likelihoods};
    std::vector<double> result(genotypes.size());
    std::transform(std::cbegin(genotypes), std::cend(genotypes), std::begin(result),
                   [&](const auto& genotype) {
                       return likelihood_model.evaluate(demote(genotype));
                   });
    maths::normalise_logs(result);
    return result;
}

LogProbabilityVector log_uniform_dist(const std::size_t n)
{
    return LogProbabilityVector(n, -std::log(static_cast<double>(n)));
}

auto generate_seeds(const std::vector<SampleName>& samples,
                    const std::vector<CancerGenotype<Haplotype>>& genotypes,
                    const LogProbabilityVector& genotype_log_priors,
                    const TumourModel::Priors& priors,
                    const HaplotypeLikelihoodCache& haplotype_log_likelihoods)
{
    std::vector<LogProbabilityVector> result{};
    result.reserve(2 + 3 * samples.size());
    result.emplace_back(genotype_log_priors);
    result.emplace_back(log_uniform_dist(genotypes.size()));
    for (const auto& sample : samples) {
        result.emplace_back(compute_germline_log_posteriors(sample, genotypes, haplotype_log_likelihoods));
        result.emplace_back(compute_germline_log_posteriors(sample, genotypes, haplotype_log_likelihoods,
                                                            priors.genotype_prior_model));
        result.emplace_back(compute_log_posteriors_with_germline_model(sample, genotypes, haplotype_log_likelihoods));
    }
    return result;
}

// Main entry point

template <std::size_t K>
TumourModel::InferredLatents
run_variational_bayes(const std::vector<SampleName>& samples,
                      std::vector<CancerGenotype<Haplotype>>&& genotypes,
                      const TumourModel::Priors& priors,
                      const HaplotypeLikelihoodCache& haplotype_log_likelihoods,
                      const VariationalBayesParameters& params)
{
    const auto prior_alphas = flatten_priors<K>(priors, samples);
    const auto genotype_log_priors = calculate_log_priors(genotypes, priors.genotype_prior_model);
    const auto log_likelihoods = compress<K>(genotypes, samples, haplotype_log_likelihoods);
    auto seeds = generate_seeds(samples, genotypes, genotype_log_priors, priors, haplotype_log_likelihoods);
    auto p = run_variational_bayes(prior_alphas, genotype_log_priors, log_likelihoods, params, seeds);
    return expand(samples, std::move(genotypes), std::move(p.first), p.second);
}

} // namespace

} // namespace model
} // namespace octopus
