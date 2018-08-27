// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "subclone_model.hpp"

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
#include "core/models/genotype/individual_model.hpp"
#include "variational_bayes_mixture_model.hpp"

namespace octopus { namespace model {

SubcloneModel::SubcloneModel(std::vector<SampleName> samples, Priors priors)
: SubcloneModel {std::move(samples), std::move(priors), AlgorithmParameters {}}
{}

SubcloneModel::SubcloneModel(std::vector<SampleName> samples, Priors priors, AlgorithmParameters parameters)
: samples_ {std::move(samples)}
, priors_ {std::move(priors)}
, parameters_ {parameters}
{}

const SubcloneModel::Priors& SubcloneModel::priors() const noexcept
{
    return priors_;
}

void SubcloneModel::prime(const std::vector<Haplotype>& haplotypes)
{
    haplotypes_ = std::addressof(haplotypes);
}

void SubcloneModel::unprime() noexcept
{
    haplotypes_ = nullptr;
}

bool SubcloneModel::is_primed() const noexcept
{
    return haplotypes_;
}

namespace {

struct IndexData
{
    const std::vector<GenotypeIndex>& genotype_indices;
    const std::vector<Haplotype>* haplotypes;
};

SubcloneModel::InferredLatents
run_variational_bayes(const std::vector<SampleName>& samples,
                      const std::vector<Genotype<Haplotype>>& genotypes,
                      const SubcloneModel::Priors& priors,
                      const HaplotypeLikelihoodCache& haplotype_log_likelihoods,
                      const SubcloneModel::AlgorithmParameters& params,
                      boost::optional<IndexData> index_data = boost::none);

} // namespace

SubcloneModel::InferredLatents
SubcloneModel::evaluate(const std::vector<Genotype<Haplotype>>& genotypes,
                   const HaplotypeLikelihoodCache& haplotype_likelihoods) const
{
    assert(!genotypes.empty());
    return run_variational_bayes(samples_, genotypes, priors_, haplotype_likelihoods, parameters_);
}

SubcloneModel::InferredLatents
SubcloneModel::evaluate(const std::vector<Genotype<Haplotype>>& genotypes,
                   const std::vector<GenotypeIndex>& genotype_indices,
                   const HaplotypeLikelihoodCache& haplotype_likelihoods) const
{
    assert(!genotypes.empty());
    assert(genotypes.size() == genotype_indices.size());
    const IndexData index_data {genotype_indices, haplotypes_};
    return run_variational_bayes(samples_, genotypes, priors_, haplotype_likelihoods, parameters_, index_data);
}

namespace {

auto evaluate_genotype_priors(const std::vector<Genotype<Haplotype>>& genotypes,
                              const SubcloneModel::Priors& priors,
                              const boost::optional<IndexData> index_data)
{
    if (index_data) {
        return evaluate(index_data->genotype_indices, priors.genotype_prior_model);
    } else {
        return evaluate(genotypes, priors.genotype_prior_model);
    }
}

template <std::size_t K>
VBAlpha<K> flatten(const SubcloneModel::Priors::GenotypeMixturesDirichletAlphas& alpha)
{
    VBAlpha<K> result {};
    std::copy_n(std::cbegin(alpha), K, std::begin(result));
    return result;
}

template <std::size_t K>
VBAlphaVector<K> flatten(const SubcloneModel::Priors::GenotypeMixturesDirichletAlphaMap& alphas,
                         const std::vector<SampleName>& samples)
{
    VBAlphaVector<K> result(samples.size());
    std::transform(std::cbegin(samples), std::cend(samples), std::begin(result),
                   [&alphas] (const auto& sample) { return flatten<K>(alphas.at(sample)); });
    return result;
}

template <std::size_t K>
VBGenotype<K>
flatten(const Genotype<Haplotype>& genotype, const SampleName& sample,
        const HaplotypeLikelihoodCache& haplotype_likelihoods)
{
    VBGenotype<K> result {};
    std::transform(std::cbegin(genotype), std::cend(genotype), std::begin(result),
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

template <std::size_t K>
SubcloneModel::Latents::GenotypeMixturesDirichletAlphas expand(VBAlpha<K>& alpha)
{
    return SubcloneModel::Latents::GenotypeMixturesDirichletAlphas(std::begin(alpha), std::end(alpha));
}

template <std::size_t K>
SubcloneModel::Latents::GenotypeMixturesDirichletAlphaMap
expand(const std::vector<SampleName>& samples, VBAlphaVector<K>&& alphas)
{
    SubcloneModel::Latents::GenotypeMixturesDirichletAlphaMap result {};
    std::transform(std::cbegin(samples), std::cend(samples), std::begin(alphas),
                   std::inserter(result, std::begin(result)),
                   [] (const auto& sample, auto&& vb_alpha) {
                       return std::make_pair(sample, expand(vb_alpha));
                   });
    return result;
}

template <std::size_t K>
SubcloneModel::InferredLatents
expand(const std::vector<SampleName>& samples, VBLatents<K>&& inferred_latents, double evidence)
{
    SubcloneModel::Latents posterior_latents {
        std::move(inferred_latents.genotype_posteriors),
        expand(samples, std::move(inferred_latents.alphas))
    };
    return {std::move(posterior_latents), evidence};
}

template <typename Container>
auto calculate_log_priors(const Container& genotypes, const GenotypePriorModel& model, const bool normalise = false)
{
    std::vector<double> result(genotypes.size());
    std::transform(std::cbegin(genotypes), std::cend(genotypes), std::begin(result),
                   [&model] (const auto& genotype) { return model.evaluate(genotype); });
    if (normalise) maths::normalise_logs(result);
    return result;
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

auto generate_exhaustive_seeds(const std::size_t n)
{
    std::vector<LogProbabilityVector> result {};
    result.reserve(n);
    for (unsigned i {0}; i < n; ++i) {
        result.push_back(make_point_seed(n, i));
    }
    return result;
}

auto generate_seeds(const std::vector<SampleName>& samples,
                    const std::vector<Genotype<Haplotype>>& genotypes,
                    const LogProbabilityVector& genotype_log_priors,
                    const HaplotypeLikelihoodCache& haplotype_log_likelihoods,
                    const SubcloneModel::Priors& priors,
                    std::size_t max_seeds,
                    boost::optional<IndexData> index_data = boost::none)
{
    if (genotypes.size() <= max_seeds) {
        return generate_exhaustive_seeds(genotypes.size());
    }
    std::vector<LogProbabilityVector> result {};
    result.reserve(1 + samples.size());
    result.push_back(genotype_log_priors);
    IndividualModel germline_model {priors.genotype_prior_model};
    for (const auto& sample : samples) {
        haplotype_log_likelihoods.prime(sample);
        IndividualModel::InferredLatents latents;
        if (index_data) {
            latents = germline_model.evaluate(genotypes, index_data->genotype_indices, haplotype_log_likelihoods);
        } else {
            latents = germline_model.evaluate(genotypes, haplotype_log_likelihoods);
        }
        result.push_back(latents.posteriors.genotype_probabilities);
        maths::log_each(result.back());
    }
    return result;
}

template <std::size_t K>
SubcloneModel::InferredLatents
run_variational_bayes(const std::vector<SampleName>& samples,
                      const std::vector<Genotype<Haplotype>>& genotypes,
                      const SubcloneModel::Priors::GenotypeMixturesDirichletAlphaMap& prior_alphas,
                      const std::vector<double>& genotype_log_priors,
                      const HaplotypeLikelihoodCache& haplotype_log_likelihoods,
                      const SubcloneModel::AlgorithmParameters& params,
                      std::vector<std::vector<double>>&& seeds)
{
    VariationalBayesParameters vb_params {params.epsilon, params.max_iterations};
    if (params.target_max_memory) {
        const auto estimated_memory_default = estimate_memory_requirement<K>(samples, haplotype_log_likelihoods, genotypes.size(), vb_params);
        if (estimated_memory_default > *params.target_max_memory) {
            vb_params.save_memory = true;
        }
    }
    const auto vb_prior_alphas = flatten<K>(prior_alphas, samples);
    const auto log_likelihoods = flatten<K>(genotypes, samples, haplotype_log_likelihoods);
    auto p = run_variational_bayes(vb_prior_alphas, genotype_log_priors, log_likelihoods, vb_params, std::move(seeds));
    return expand(samples, std::move(p.first), p.second);
}

SubcloneModel::InferredLatents
run_variational_bayes_helper(const std::vector<SampleName>& samples,
                             const std::vector<Genotype<Haplotype>>& genotypes,
                             const SubcloneModel::Priors::GenotypeMixturesDirichletAlphaMap& prior_alphas,
                             std::vector<double> genotype_log_priors,
                             const HaplotypeLikelihoodCache& haplotype_log_likelihoods,
                             const SubcloneModel::AlgorithmParameters& params,
                             std::vector<std::vector<double>>&& seeds)
{
    using std::move;
    switch (genotypes.front().ploidy()) {
    case 2: return run_variational_bayes<2>(samples, genotypes, prior_alphas, move(genotype_log_priors),
                                            haplotype_log_likelihoods, params, move(seeds));
    case 3: return run_variational_bayes<3>(samples, genotypes, prior_alphas, move(genotype_log_priors),
                                            haplotype_log_likelihoods, params, move(seeds));
    case 4: return run_variational_bayes<4>(samples, genotypes, prior_alphas, move(genotype_log_priors),
                                            haplotype_log_likelihoods, params, move(seeds));
    case 5: return run_variational_bayes<5>(samples, genotypes, prior_alphas, move(genotype_log_priors),
                                            haplotype_log_likelihoods, params, move(seeds));
    case 6: return run_variational_bayes<6>(samples, genotypes, prior_alphas, move(genotype_log_priors),
                                            haplotype_log_likelihoods, params, move(seeds));
    case 7: return run_variational_bayes<7>(samples, genotypes, prior_alphas, move(genotype_log_priors),
                                            haplotype_log_likelihoods, params, move(seeds));
    case 8: return run_variational_bayes<8>(samples, genotypes, prior_alphas, move(genotype_log_priors),
                                            haplotype_log_likelihoods, params, move(seeds));
    case 9: return run_variational_bayes<9>(samples, genotypes, prior_alphas, move(genotype_log_priors),
                                            haplotype_log_likelihoods, params, move(seeds));
    case 10: return run_variational_bayes<10>(samples, genotypes, prior_alphas, move(genotype_log_priors),
                                              haplotype_log_likelihoods, params, move(seeds));
    default: throw UnimplementedFeatureError {"ploidies above 10", "SubcloneModel"};
    }
}

// Main entry point

SubcloneModel::InferredLatents
run_variational_bayes(const std::vector<SampleName>& samples,
                      const std::vector<Genotype<Haplotype>>& genotypes,
                      const SubcloneModel::Priors& priors,
                      const HaplotypeLikelihoodCache& haplotype_log_likelihoods,
                      const SubcloneModel::AlgorithmParameters& params,
                      const boost::optional<IndexData> index_data)
{
    auto genotype_log_priors = evaluate_genotype_priors(genotypes, priors, index_data);
    auto seeds = generate_seeds(samples, genotypes, genotype_log_priors, haplotype_log_likelihoods, priors, params.max_seeds, index_data);
    return run_variational_bayes_helper(samples, genotypes, priors.alphas, std::move(genotype_log_priors),
                                        haplotype_log_likelihoods, params, std::move(seeds));
}

} // namespace

} // namespace model
} // namespace octopus
