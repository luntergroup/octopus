// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "cnv_model.hpp"

#include <algorithm>
#include <numeric>
#include <iterator>
#include <functional>
#include <cstddef>
#include <cmath>
#include <cassert>

#include "utils/maths.hpp"
#include "logging/logging.hpp"
#include "core/models/genotype/individual_model.hpp"
#include "variational_bayes_mixture_model.hpp"

namespace octopus { namespace model {

CNVModel::CNVModel(std::vector<SampleName> samples, const unsigned ploidy, Priors priors)
: CNVModel {std::move(samples), ploidy, std::move(priors), AlgorithmParameters {}}
{}

CNVModel::CNVModel(std::vector<SampleName> samples, const unsigned ploidy, Priors priors,
                   AlgorithmParameters parameters)
: samples_ {std::move(samples)}
, ploidy_ {ploidy}
, priors_ {std::move(priors)}
, parameters_ {parameters}
{}

const CNVModel::Priors& CNVModel::priors() const noexcept
{
    return priors_;
}

namespace {

template <std::size_t K>
CNVModel::InferredLatents
run_variational_bayes(const std::vector<SampleName>& samples,
                      std::vector<Genotype<Haplotype>>&& genotypes,
                      const CNVModel::Priors& priors,
                      const HaplotypeLikelihoodCache& haplotype_log_likelihoods,
                      const VariationalBayesParameters& params);
    
} // namespace

CNVModel::InferredLatents
CNVModel::evaluate(std::vector<Genotype<Haplotype>> genotypes,
                   const HaplotypeLikelihoodCache& haplotype_likelihoods) const
{
    assert(!genotypes.empty());
    assert(ploidy_ < 4);
    const VariationalBayesParameters vb_params {parameters_.epsilon, parameters_.max_iterations};
    switch (ploidy_) {
        case 1: return run_variational_bayes<1>(samples_, std::move(genotypes), priors_,
                                                haplotype_likelihoods, vb_params);
        case 2: return run_variational_bayes<2>(samples_, std::move(genotypes), priors_,
                                                haplotype_likelihoods, vb_params);
        default: return run_variational_bayes<3>(samples_, std::move(genotypes), priors_,
                                                 haplotype_likelihoods, vb_params);
    }
}

namespace {

template <std::size_t K>
VBAlpha<K> flatten(const CNVModel::Priors::GenotypeMixturesDirichletAlphas& alpha)
{
    VBAlpha<K> result;
    std::copy_n(std::cbegin(alpha), K, std::begin(result));
    return result;
}

template <std::size_t K>
VBAlphaVector<K> flatten(const CNVModel::Priors& priors, const std::vector<SampleName>& samples)
{
    VBAlphaVector<K> result(samples.size());
    std::transform(std::cbegin(samples), std::cend(samples), std::begin(result),
                   [&priors] (const auto& sample) {
                       return flatten<K>(priors.alphas.at(sample));
                   });
    return result;
}

template <std::size_t K>
VBGenotype<K>
flatten(const Genotype<Haplotype>& genotype, const SampleName& sample,
        const HaplotypeLikelihoodCache& haplotype_likelihoods)
{
    VBGenotype<K> result;
    std::transform(std::cbegin(genotype), std::cend(genotype),
                   std::begin(result),
                   [&sample, &haplotype_likelihoods] (const Haplotype& haplotype)
                   -> std::reference_wrapper<const VBReadLikelihoodArray::BaseType> {
                       return std::cref(haplotype_likelihoods(sample, haplotype));
                   });
    return result;
}

template <std::size_t K>
VBGenotypeVector<K>
flatten(const std::vector<Genotype<Haplotype>>& genotypes, const SampleName& sample,
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
flatten(const std::vector<Genotype<Haplotype>>& genotypes,
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

template <typename Container>
std::vector<double> calculate_log_priors(const Container& genotypes, const GenotypePriorModel& model)
{
    std::vector<double> result(genotypes.size());
    std::transform(std::cbegin(genotypes), std::cend(genotypes), std::begin(result),
                   [&model](const auto& genotype) {
                       return model.evaluate(genotype);
                   });
    maths::normalise_logs(result);
    return result;
}

CNVModel::Latents::GenotypeProbabilityMap
expand(std::vector<Genotype<Haplotype>>&& genotypes, LogProbabilityVector&& genotype_log_posteriors)
{
    CNVModel::Latents::GenotypeProbabilityMap result {};
    std::transform(std::make_move_iterator(std::begin(genotypes)),
                   std::make_move_iterator(std::end(genotypes)),
                   std::begin(genotype_log_posteriors),
                   std::inserter(result, std::begin(result)),
                   [] (auto&& g, auto p) {
                       return std::make_pair(std::move(g), p);
                   });
    return result;
}

template <std::size_t K>
CNVModel::Latents::GenotypeMixturesDirichletAlphas expand(VBAlpha<K>& alpha)
{
    return CNVModel::Latents::GenotypeMixturesDirichletAlphas(std::begin(alpha), std::end(alpha));
}

template <std::size_t K>
CNVModel::Latents::GenotypeMixturesDirichletAlphaMap
expand(const std::vector<SampleName>& samples, VBAlphaVector<K>&& alphas)
{
    CNVModel::Latents::GenotypeMixturesDirichletAlphaMap result {};
    std::transform(std::cbegin(samples), std::cend(samples), std::begin(alphas),
                   std::inserter(result, std::begin(result)),
                   [] (const auto& sample, auto&& vb_alpha) {
                       return std::make_pair(sample, expand(vb_alpha));
                   });
    return result;
}

template <std::size_t K>
CNVModel::InferredLatents
expand(const std::vector<SampleName>& samples, std::vector<Genotype<Haplotype>>&& genotypes,
       VBLatents<K>&& inferred_latents, double evidence)
{
    CNVModel::Latents posterior_latents {
        expand(std::move(genotypes), std::move(inferred_latents.genotype_posteriors)),
        expand(samples, std::move(inferred_latents.alphas))
    };
    return CNVModel::InferredLatents {std::move(posterior_latents), evidence};
}

LogProbabilityVector log_uniform_dist(const std::size_t n)
{
    return LogProbabilityVector(n, -std::log(static_cast<double>(n)));
}

auto generate_seeds(const std::vector<SampleName>& samples,
                    const std::vector<Genotype<Haplotype>>& genotypes,
                    const LogProbabilityVector& genotype_log_priors,
                    const CNVModel::Priors& priors,
                    const HaplotypeLikelihoodCache& haplotype_log_likelihoods)
{
    std::vector<LogProbabilityVector> result {};
    result.reserve(2 + 2 * samples.size());
    result.push_back(genotype_log_priors);
    result.push_back(log_uniform_dist(genotypes.size()));
    model::IndividualModel germline_model {priors.genotype_prior_model};
    for (const auto& sample : samples) {
        haplotype_log_likelihoods.prime(sample);
        const auto latents = germline_model.evaluate(genotypes, haplotype_log_likelihoods);
        result.push_back(latents.posteriors.genotype_probabilities);
        maths::log_each(result.back());
        result.push_back(latents.posteriors.genotype_probabilities);
        for (auto& p : result.back()) p = 1.0 - p;
        maths::log_each(result.back());
    }
    return result;
}

// Main entry point

template <std::size_t K>
CNVModel::InferredLatents
run_variational_bayes(const std::vector<SampleName>& samples,
                      std::vector<Genotype<Haplotype>>&& genotypes,
                      const CNVModel::Priors& priors,
                      const HaplotypeLikelihoodCache& haplotype_log_likelihoods,
                      const VariationalBayesParameters& params)
{
    const auto prior_alphas = flatten<K>(priors, samples);
    const auto genotype_log_priors = calculate_log_priors(genotypes, priors.genotype_prior_model);
    const auto log_likelihoods = flatten<K>(genotypes, samples, haplotype_log_likelihoods);
    auto seeds = generate_seeds(samples, genotypes, genotype_log_priors,
                                priors, haplotype_log_likelihoods);
    auto p = run_variational_bayes(prior_alphas, genotype_log_priors,
                                   log_likelihoods, params, std::move(seeds));
    return expand(samples, std::move(genotypes), std::move(p.first), p.second);
}

} // namespace

} // namespace model
} // namespace octopus
