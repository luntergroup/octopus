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

#include "exceptions/unimplemented_feature_error.hpp"
#include "utils/maths.hpp"
#include "logging/logging.hpp"
#include "germline_likelihood_model.hpp"
#include "fixed_mixture_genotype_likelihood_model.hpp"
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

TumourModel::InferredLatents
run_variational_bayes(const std::vector<SampleName>& samples,
                      const std::vector<CancerGenotype<Haplotype>>& genotypes,
                      const TumourModel::Priors& priors,
                      const HaplotypeLikelihoodCache& haplotype_log_likelihoods,
                      const TumourModel::AlgorithmParameters& params);

TumourModel::InferredLatents
run_variational_bayes(const std::vector<SampleName>& samples,
                      const std::vector<CancerGenotype<Haplotype>>& genotypes,
                      const std::vector<CancerGenotypeIndex>& genotype_indices,
                      const TumourModel::Priors& priors,
                      const HaplotypeLikelihoodCache& haplotype_log_likelihoods,
                      const TumourModel::AlgorithmParameters& params);

} // namespace

TumourModel::InferredLatents
TumourModel::evaluate(const std::vector<CancerGenotype<Haplotype>>& genotypes,
                      const HaplotypeLikelihoodCache& haplotype_likelihoods) const
{
    assert(!genotypes.empty());
    return run_variational_bayes(samples_, genotypes, priors_, haplotype_likelihoods, parameters_);
}

TumourModel::InferredLatents
TumourModel::evaluate(const std::vector<CancerGenotype<Haplotype>>& genotypes,
                      const std::vector<CancerGenotypeIndex>& genotype_indices,
                      const HaplotypeLikelihoodCache& haplotype_likelihoods) const
{
    assert(!genotypes.empty());
    assert(genotypes.size() == genotype_indices.size());
    return run_variational_bayes(samples_, genotypes, genotype_indices, priors_, haplotype_likelihoods, parameters_);
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
auto copy_cref(const Genotype<Haplotype>& genotype, const SampleName& sample,
               const HaplotypeLikelihoodCache& haplotype_likelihoods,
               typename VBGenotype<K>::iterator result_itr)
{
    return std::transform(std::cbegin(genotype), std::cend(genotype), result_itr,
                          [&sample, &haplotype_likelihoods] (const Haplotype& haplotype)
                          -> std::reference_wrapper<const VBReadLikelihoodArray::BaseType> {
                            return std::cref(haplotype_likelihoods(sample, haplotype));
                           });
}

template <std::size_t K>
VBGenotype<K>
flatten(const CancerGenotype<Haplotype>& genotype, const SampleName& sample,
        const HaplotypeLikelihoodCache& haplotype_likelihoods)
{
    VBGenotype<K> result {};
    assert(genotype.ploidy() == K);
    auto itr = copy_cref<K>(genotype.germline(), sample, haplotype_likelihoods, std::begin(result));
    copy_cref<K>(genotype.somatic(), sample, haplotype_likelihoods, itr);
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
expand(const std::vector<SampleName>& samples, VBLatents<K>&& inferred_latents,
       std::vector<double> genotype_log_priors, double evidence)
{
    TumourModel::Latents posterior_latents {std::move(inferred_latents.genotype_posteriors),
                                            expand(samples, std::move(inferred_latents.alphas))};
    return {std::move(posterior_latents), std::move(genotype_log_priors), evidence};
}

LogProbabilityVector log_uniform_dist(const std::size_t n)
{
    return LogProbabilityVector(n, -std::log(static_cast<double>(n)));
}

auto make_point_seed(const std::size_t num_genotypes, const std::size_t idx, const double p = 0.9999)
{
    LogProbabilityVector result(num_genotypes, num_genotypes > 1 ? std::log((1 - p) / (num_genotypes - 1)) : 0);
    if (num_genotypes > 1) result[idx] = std::log(p);
    return result;
}

void make_point_seeds(const std::size_t num_genotypes, const std::vector<std::size_t>& indices,
                      std::vector<LogProbabilityVector>& result, const double p = 0.9999)
{
    result.reserve(result.size() + indices.size());
    std::transform(std::cbegin(indices), std::cend(indices), std::back_inserter(result),
                   [=] (auto idx) { return make_point_seed(num_genotypes, idx, p); });
}

auto make_multipoint_seed(const std::size_t num_genotypes, const std::vector<std::size_t>& indices,
                          double p = 0.9999999)
{
    assert(num_genotypes >= indices.size());
    LogProbabilityVector result(num_genotypes, num_genotypes > 1 ? std::log((1 - p) / (num_genotypes - 1)) : 0);
    if (num_genotypes > 1) {
        p /= indices.size();
        for (auto idx : indices) result[idx] = std::log(p);
    }
    return result;
}

auto make_range_seed(const std::size_t num_genotypes, const std::size_t begin, const std::size_t n,
                     const double p = 0.9999999)
{
    LogProbabilityVector result(num_genotypes, std::log((1 - p) / (num_genotypes - n)));
    std::fill_n(std::next(std::begin(result), begin), n, std::log(p / n));
    return result;
}

auto make_range_seed(const std::vector<CancerGenotype<Haplotype>>& genotypes, const Genotype<Haplotype>& germline,
                     const double p = 0.9999999)
{
    auto itr1 = std::find_if(std::cbegin(genotypes), std::cend(genotypes), [&] (const auto& g) { return g.germline() == germline; });
    auto itr2 = std::find_if_not(std::next(itr1), std::cend(genotypes), [&] (const auto& g) { return g.germline() == germline; });
    return make_range_seed(genotypes.size(), std::distance(std::cbegin(genotypes), itr1), std::distance(itr1, itr2));
}

namespace debug {

template <typename S>
void print_top(S&& stream, const std::vector<CancerGenotype<Haplotype>>& genotypes,
               const LogProbabilityVector& probs, std::size_t n)
{
    assert(probs.size() == genotypes.size());
    n = std::min(n, genotypes.size());
    std::vector<std::pair<CancerGenotype<Haplotype>, double> > pairs {};
    pairs.reserve(genotypes.size());
    std::transform(std::cbegin(genotypes), std::cend(genotypes), std::cbegin(probs), std::back_inserter(pairs),
                   [] (const auto& g, auto p) { return std::make_pair(g, p); });
    const auto mth = std::next(std::begin(pairs), n);
    std::partial_sort(std::begin(pairs), mth, std::end(pairs),
                      [] (const auto& lhs, const auto& rhs) { return lhs.second > rhs.second; });
    std::for_each(std::begin(pairs), mth, [&] (const auto& p) {
        octopus::debug::print_variant_alleles(stream, p.first);
        stream << " " << p.second << '\n';
    });
}

void print_top(const std::vector<CancerGenotype<Haplotype>>& genotypes,
               const LogProbabilityVector& probs, std::size_t n = 10)
{
    print_top(std::cout, genotypes, probs, n);
}

} // namespace debug

bool is_somatic_expected(const SampleName& sample, const TumourModel::Priors& priors)
{
    const auto& alphas = priors.alphas.at(sample);
    auto e = maths::dirichlet_expectation(alphas.size() - 1, alphas);
    return e > 0.05;
}

void add_to(const LogProbabilityVector& other, LogProbabilityVector& result)
{
    std::transform(std::cbegin(other), std::cend(other), std::cbegin(result), std::begin(result),
                   [] (auto a, auto b) { return a + b; });
}

auto generate_exhaustive_seeds(const std::size_t n)
{
    std::vector<LogProbabilityVector> result {};
    result.reserve(n);
    for (unsigned i {0}; i < n; ++i) {
        result.push_back(make_point_seed(n, i));
    }
    return result;
}

std::vector<LogProbabilityVector>
compute_genotype_likelihoods_with_fixed_mixture_model(const std::vector<SampleName>& samples,
                                                      const std::vector<CancerGenotype<Haplotype>>& genotypes,
                                                      const HaplotypeLikelihoodCache& haplotype_log_likelihoods,
                                                      const TumourModel::Priors::GenotypeMixturesDirichletAlphaMap& priors)
{
    FixedMixtureGenotypeLikelihoodModel model {haplotype_log_likelihoods};
    std::vector<LogProbabilityVector> result {};
    result.reserve(samples.size());
    for (const auto& sample : samples) {
        const auto& sample_priors = priors.at(sample);
        model.set_mixtures(maths::dirichlet_expectation(sample_priors));
        model.cache().prime(sample);
        result.push_back(evaluate(genotypes, model));
    }
    return result;
}

double evaluate(const CancerGenotype<Haplotype>& genotype, const GermlineLikelihoodModel& model)
{
    return model.evaluate(demote(genotype));
}

LogProbabilityVector evaluate(const std::vector<CancerGenotype<Haplotype>>& genotypes, const GermlineLikelihoodModel& model)
{
    LogProbabilityVector result(genotypes.size());
    std::transform(std::cbegin(genotypes), std::cend(genotypes), std::begin(result),
                   [&] (const auto& genotype) { return evaluate(genotype, model); });
    return result;
}

std::vector<LogProbabilityVector>
compute_genotype_likelihoods_with_germline_model(const std::vector<SampleName>& samples,
                                                 const std::vector<CancerGenotype<Haplotype>>& genotypes,
                                                 const HaplotypeLikelihoodCache& haplotype_log_likelihoods)
{
    const GermlineLikelihoodModel model {haplotype_log_likelihoods};
    std::vector<LogProbabilityVector> result {};
    result.reserve(samples.size());
    for (const auto& sample : samples) {
        model.cache().prime(sample);
        result.push_back(evaluate(genotypes, model));
    }
    return result;
}

LogProbabilityVector evaluate_germlines(const std::vector<CancerGenotype<Haplotype>>& genotypes, const GermlineLikelihoodModel& model)
{
    std::unordered_map<Genotype<Haplotype>, double> cache {};
    cache.reserve(genotypes.size());
    LogProbabilityVector result(genotypes.size());
    std::transform(std::cbegin(genotypes), std::cend(genotypes), std::begin(result),
                   [&] (const auto& genotype) {
        const auto cache_itr = cache.find(genotype.germline());
        if (cache_itr != std::cend(cache)) return cache_itr->second;
        const auto result = model.evaluate(genotype.germline());
        cache.emplace(genotype.germline(), result);
        return result;
    });
    return result;
}

std::vector<LogProbabilityVector>
compute_germline_genotype_likelihoods_with_germline_model(const std::vector<SampleName>& samples,
                                                          const std::vector<CancerGenotype<Haplotype>>& genotypes,
                                                          const HaplotypeLikelihoodCache& haplotype_log_likelihoods)
{
    const GermlineLikelihoodModel model {haplotype_log_likelihoods};
    std::vector<LogProbabilityVector> result {};
    result.reserve(samples.size());
    for (const auto& sample : samples) {
        model.cache().prime(sample);
        result.push_back(evaluate_germlines(genotypes, model));
    }
    return result;
}

auto add_all(const std::vector<ProbabilityVector>& likelihoods)
{
    assert(!likelihoods.empty());
    ProbabilityVector result(likelihoods.front().size());
    for (const auto& probabilities : likelihoods) {
        add_to(probabilities, result);
    }
    return result;
}

auto add_all_and_normalise(const std::vector<LogProbabilityVector>& log_likelihoods)
{
    auto result = add_all(log_likelihoods);
    maths::normalise_logs(result);
    return result;
}

auto add(const LogProbabilityVector& lhs, const LogProbabilityVector& rhs)
{
    auto result = lhs;
    add_to(rhs, result);
    return result;
}

auto add_and_normalise(const LogProbabilityVector& lhs, const LogProbabilityVector& rhs)
{
    auto result = add(lhs, rhs);
    maths::normalise_logs(result);
    return result;
}

auto generate_seeds(const std::vector<SampleName>& samples,
                    const std::vector<CancerGenotype<Haplotype>>& genotypes,
                    const LogProbabilityVector& genotype_log_priors,
                    const HaplotypeLikelihoodCache& haplotype_log_likelihoods,
                    const TumourModel::Priors& priors,
                    std::size_t max_seeds)
{
    if (genotypes.size() <= max_seeds) {
        return generate_exhaustive_seeds(genotypes.size());
    }
    std::vector<LogProbabilityVector> result {};
    result.reserve(max_seeds);
    if (max_seeds == 0) return result;
    auto sample_prior_mixture_likelihoods = compute_genotype_likelihoods_with_fixed_mixture_model(samples, genotypes, haplotype_log_likelihoods, priors.alphas);
    auto prior_mixture_likelihoods = add_all_and_normalise(sample_prior_mixture_likelihoods);
    auto prior_mixture_posteriors = add_and_normalise(genotype_log_priors, prior_mixture_likelihoods);
    result.push_back(prior_mixture_posteriors); // 1
    --max_seeds;
    if (max_seeds == 0) return result;
    auto sample_normal_likelihoods = compute_genotype_likelihoods_with_germline_model(samples, genotypes, haplotype_log_likelihoods);
    auto normal_likelihoods = add_all_and_normalise(sample_prior_mixture_likelihoods);
    auto normal_posteriors = add_and_normalise(genotype_log_priors, normal_likelihoods);
    result.push_back(normal_posteriors); // 2
    --max_seeds;
    if (max_seeds == 0) return result;
    result.push_back(prior_mixture_likelihoods); // 3
    if (max_seeds == 0) return result;
    result.push_back(normal_likelihoods); // 4
    if (max_seeds == 0) return result;
    auto combined_model_likelihoods = add_and_normalise(prior_mixture_likelihoods, normal_likelihoods);
    auto combined_model_posteriors = add_and_normalise(genotype_log_priors, combined_model_likelihoods);
    result.push_back(combined_model_posteriors); // 5
    if (max_seeds == 0) return result;
    result.push_back(combined_model_likelihoods); // 6
    if (max_seeds == 0) return result;
    auto sample_germline_likelihoods = compute_germline_genotype_likelihoods_with_germline_model(samples, genotypes, haplotype_log_likelihoods);
    result.push_back(add_all_and_normalise(sample_germline_likelihoods)); // 7
    --max_seeds;
    if (max_seeds == 0) return result;
//    if (samples.size() > 1) {
//        for (auto& log_likelihoods : sample_prior_mixture_likelihoods) {
//            result.push_back(add_and_normalise(genotype_log_priors, log_likelihoods));
//            result.push_back(log_likelihoods);
//            if (max_seeds == 0) return result;
//            maths::normalise_logs(log_likelihoods);
//            result.push_back(std::move(log_likelihoods));
//            if (max_seeds == 0) return result;
//        }
//        for (auto& log_likelihoods : sample_normal_likelihoods) {
//            result.push_back(add_and_normalise(genotype_log_priors, log_likelihoods));
//            result.push_back(log_likelihoods);
//            if (max_seeds == 0) return result;
//            maths::normalise_logs(log_likelihoods);
//            result.push_back(std::move(log_likelihoods));
//            if (max_seeds == 0) return result;
//        }
//    }
    result.push_back(genotype_log_priors); // 8
    maths::normalise_logs(result.back());
    if (max_seeds == 0) return result;
    std::vector<std::pair<double, std::size_t>> ranked_prior_mixture_posteriors {};
    ranked_prior_mixture_posteriors.reserve(genotypes.size());
    for (std::size_t i {0}; i < genotypes.size(); ++i) {
        ranked_prior_mixture_posteriors.push_back({prior_mixture_posteriors[i], i});
    }
    const auto nth_ranked_prior_mixture_posteriors_itr = std::next(std::begin(ranked_prior_mixture_posteriors), max_seeds);
    std::partial_sort(std::begin(ranked_prior_mixture_posteriors),
                      nth_ranked_prior_mixture_posteriors_itr,
                      std::end(ranked_prior_mixture_posteriors), std::greater<> {});
    std::transform(std::begin(ranked_prior_mixture_posteriors), nth_ranked_prior_mixture_posteriors_itr,
                   std::back_inserter(result), [&] (const auto& p) { return make_point_seed(genotypes.size(), p.second); });
    return result;
}

template <std::size_t K>
TumourModel::InferredLatents
run_variational_bayes(const std::vector<SampleName>& samples,
                      const std::vector<CancerGenotype<Haplotype>>& genotypes,
                      const TumourModel::Priors::GenotypeMixturesDirichletAlphaMap& prior_alphas,
                      std::vector<double> genotype_log_priors,
                      const HaplotypeLikelihoodCache& haplotype_log_likelihoods,
                      const VariationalBayesParameters& params,
                      std::vector<std::vector<double>>&& seeds)
{
    const auto vb_prior_alphas = flatten<K>(prior_alphas, samples);
    const auto log_likelihoods = flatten<K>(genotypes, samples, haplotype_log_likelihoods);
    auto p = run_variational_bayes(vb_prior_alphas, genotype_log_priors, log_likelihoods, params, std::move(seeds));
    return expand(samples, std::move(p.first), std::move(genotype_log_priors), p.second);
}

TumourModel::InferredLatents
run_variational_bayes_helper(const std::vector<SampleName>& samples,
                             const std::vector<CancerGenotype<Haplotype>>& genotypes,
                             const TumourModel::Priors::GenotypeMixturesDirichletAlphaMap& prior_alphas,
                             std::vector<double> genotype_log_priors,
                             const HaplotypeLikelihoodCache& haplotype_log_likelihoods,
                             const TumourModel::AlgorithmParameters& params,
                             std::vector<std::vector<double>>&& seeds)
{
    const VariationalBayesParameters vb_params {params.epsilon, params.max_iterations};
    using std::move;
    switch (genotypes.front().ploidy()) {
        case 2: return run_variational_bayes<2>(samples, genotypes, prior_alphas, move(genotype_log_priors),
                                                haplotype_log_likelihoods, vb_params, move(seeds));
        case 3: return run_variational_bayes<3>(samples, genotypes, prior_alphas, move(genotype_log_priors),
                                                haplotype_log_likelihoods, vb_params, move(seeds));
        case 4: return run_variational_bayes<4>(samples, genotypes, prior_alphas, move(genotype_log_priors),
                                                haplotype_log_likelihoods, vb_params, move(seeds));
        case 5: return run_variational_bayes<5>(samples, genotypes, prior_alphas, move(genotype_log_priors),
                                                haplotype_log_likelihoods, vb_params, move(seeds));
        case 6: return run_variational_bayes<6>(samples, genotypes, prior_alphas, move(genotype_log_priors),
                                                haplotype_log_likelihoods, vb_params, move(seeds));
        case 7: return run_variational_bayes<7>(samples, genotypes, prior_alphas, move(genotype_log_priors),
                                                haplotype_log_likelihoods, vb_params, move(seeds));
        case 8: return run_variational_bayes<8>(samples, genotypes, prior_alphas, move(genotype_log_priors),
                                                haplotype_log_likelihoods, vb_params, move(seeds));
        default: throw UnimplementedFeatureError {"tumour model ploidies above 8", "TumourModel"};
    }
}

// Main entry point

TumourModel::InferredLatents
run_variational_bayes(const std::vector<SampleName>& samples,
                      const std::vector<CancerGenotype<Haplotype>>& genotypes,
                      const TumourModel::Priors& priors,
                      const HaplotypeLikelihoodCache& haplotype_log_likelihoods,
                      const TumourModel::AlgorithmParameters& params)
{
    auto genotype_log_priors = evaluate(genotypes, priors.genotype_prior_model);
    auto seeds = generate_seeds(samples, genotypes, genotype_log_priors, haplotype_log_likelihoods, priors, params.max_seeds);
    return run_variational_bayes_helper(samples, genotypes, priors.alphas, std::move(genotype_log_priors),
                                        haplotype_log_likelihoods, params, std::move(seeds));
}

TumourModel::InferredLatents
run_variational_bayes(const std::vector<SampleName>& samples,
                      const std::vector<CancerGenotype<Haplotype>>& genotypes,
                      const std::vector<CancerGenotypeIndex>& genotype_indices,
                      const TumourModel::Priors& priors,
                      const HaplotypeLikelihoodCache& haplotype_log_likelihoods,
                      const TumourModel::AlgorithmParameters& params)
{
    auto genotype_log_priors = evaluate(genotype_indices, priors.genotype_prior_model);
    auto seeds = generate_seeds(samples, genotypes, genotype_log_priors, haplotype_log_likelihoods, priors, params.max_seeds);
    return run_variational_bayes_helper(samples, genotypes, priors.alphas, std::move(genotype_log_priors),
                                        haplotype_log_likelihoods, params, std::move(seeds));
}

} // namespace

} // namespace model
} // namespace octopus
