// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "single_cell_model.hpp"

#include <utility>
#include <cassert>
#include <iterator>
#include <algorithm>
#include <numeric>

#include "utils/select_top_k.hpp"
#include "utils/maths.hpp"
#include "subclone_model.hpp"
#include "population_model.hpp"
#include "uniform_population_prior_model.hpp"
#include "individual_model.hpp"

namespace octopus { namespace model {

SingleCellModel::SingleCellModel(std::vector<SampleName> samples, SingleCellPriorModel prior_model,
                                 Parameters parameters, AlgorithmParameters config)
: samples_ {std::move(samples)}
, prior_model_ {std::move(prior_model)}
, parameters_ {std::move(parameters)}
, config_ {std::move(config)}
{}

SingleCellModel::Inferences
SingleCellModel::evaluate(const std::vector<Genotype<Haplotype>>& genotypes,
                          const HaplotypeLikelihoodArray& haplotype_likelihoods) const
{
    Inferences result {};
    if (prior_model_.phylogeny().size() == 1) {
        SubcloneModel::Priors subclone_priors {prior_model_.germline_prior_model(), {}};
        const auto ploidy = genotypes.front().ploidy();
        for (const auto& sample : samples_) {
            subclone_priors.alphas.emplace(sample, SubcloneModel::Priors::GenotypeMixturesDirichletAlphas(ploidy, parameters_.dropout_concentration));
        }
        SubcloneModel helper_model {samples_, std::move(subclone_priors)};
        auto subclone_inferences = helper_model.evaluate(genotypes, haplotype_likelihoods);
        Inferences::GroupInferences founder {};
        founder.genotype_posteriors = std::move(subclone_inferences.posteriors.genotype_probabilities);
        founder.sample_attachment_posteriors.assign(samples_.size(), 1.0);
        result.phylogeny.set_founder({0, std::move(founder)});
        result.log_evidence = subclone_inferences.approx_log_evidence;
    } else {
        const auto genotype_combinations = propose_genotype_combinations(genotypes, haplotype_likelihoods);
        const auto genotype_combination_priors = calculate_genotype_priors(genotype_combinations, genotypes);
        const auto vb_haplotype_likelihoods = make_likelihood_matrix(genotype_combinations, genotypes, haplotype_likelihoods);
        auto seeds = propose_seeds(genotype_combinations, genotypes, genotype_combination_priors, haplotype_likelihoods);
        auto vb_inferences = posterior_model_.evaluate(genotype_combination_priors, vb_haplotype_likelihoods,
                                                       parameters_.group_concentration, parameters_.dropout_concentration,
                                                       std::move(seeds));
        for (std::size_t group_idx {0}; group_idx < prior_model_.phylogeny().size(); ++group_idx) {
            Inferences::GroupInferences group {};
            group.sample_attachment_posteriors.resize(samples_.size());
            for (std::size_t sample_idx {0}; sample_idx < samples_.size(); ++sample_idx) {
                group.sample_attachment_posteriors[sample_idx] = vb_inferences.group_responsabilities[sample_idx][group_idx];
            }
            // Marginalise over genotypes
            group.genotype_posteriors.resize(genotypes.size());
            for (std::size_t genotype_combo_idx {0}; genotype_combo_idx < genotype_combinations.size(); ++genotype_combo_idx) {
                group.genotype_posteriors[genotype_combinations[genotype_combo_idx][group_idx]] += vb_inferences.genotype_posteriors[genotype_combo_idx];
            }
            if (group_idx == 0) {
                result.phylogeny.set_founder({group_idx, std::move(group)});
            } else {
                const auto ancestor_idx = prior_model_.phylogeny().ancestor(group_idx).id;
                result.phylogeny.add_descendant({group_idx, std::move(group)}, ancestor_idx);
            }
        }
        result.log_evidence = vb_inferences.approx_log_evidence;
    }
    return result;
}

SingleCellModel::Inferences
SingleCellModel::evaluate(const std::vector<GenotypeIndex>& genotypes,
                          const HaplotypeLikelihoodArray& haplotype_likelihoods) const
{
    return {};
}

// private methods

namespace {

auto log(std::size_t base, std::size_t x)
{
    return std::log2(x) / std::log2(base);
}

auto num_combinations(const std::size_t num_genotypes, const std::size_t num_samples)
{
    static constexpr auto max_combinations = std::numeric_limits<std::size_t>::max();
    if (num_samples <= log(num_genotypes, max_combinations)) {
        return static_cast<std::size_t>(std::pow(num_genotypes, num_samples));
    } else {
        return max_combinations;
    }
}

auto pool_likelihood(const std::vector<SampleName>& samples,
                     const std::vector<Haplotype>& haplotypes,
                     const HaplotypeLikelihoodArray& haplotype_likelihoods)
{
    static const SampleName pooled_sample {"pool"};
    auto result = merge_samples(samples, pooled_sample, haplotypes, haplotype_likelihoods);
    result.prime(pooled_sample);
    return result;
}

auto kl_divergence(const std::vector<double>& p, const std::vector<double>& q) noexcept
{
    return std::inner_product(std::cbegin(p), std::cend(p), std::cbegin(q), 0.0,
                              std::plus<> {}, [] (const auto a, const auto b) {
        return (a > 0 && b > 0) ? a * std::log(a / b) : 0.0;
    });
}

auto symmetric_kl_divergence(const std::vector<double>& p, const std::vector<double>& q) noexcept
{
    return kl_divergence(p, q) + kl_divergence(q, p);
}

} // namespace

SingleCellModel::GenotypeCombinationVector
SingleCellModel::propose_genotype_combinations(const std::vector<Genotype<Haplotype>>& genotypes,
                                               const HaplotypeLikelihoodArray& haplotype_likelihoods) const
{
    const auto num_groups = prior_model_.phylogeny().size();
    const auto max_possible_combinations = num_combinations(genotypes.size(), num_groups);
    if (max_possible_combinations <= config_.max_genotype_combinations) {
        return propose_all_genotype_combinations(genotypes);
    } else {
        // 1. Run population model
        // 2. Cluster samples
        // 3. Run individual model on merged reads
        // 4. Select top combinations using cluster marginal posteriors
        
        UniformPopulationPriorModel population_prior_model {};
        PopulationModel::Options population_model_options {};
        population_model_options.max_joint_genotypes = config_.max_genotype_combinations;
        PopulationModel population_model {population_prior_model, population_model_options};
        const auto population_inferences = population_model.evaluate(samples_, genotypes, haplotype_likelihoods);
        const auto& population_genotype_posteriors = population_inferences.posteriors.marginal_genotype_probabilities;
        
        std::vector<std::vector<double>> kl_divergences(samples_.size(), std::vector<double>(samples_.size()));
        for (std::size_t lhs_sample_idx {0}; lhs_sample_idx < samples_.size(); ++lhs_sample_idx) {
            for (std::size_t rhs_sample_idx {lhs_sample_idx + 1}; rhs_sample_idx < samples_.size(); ++rhs_sample_idx) {
                auto kl = symmetric_kl_divergence(population_genotype_posteriors[lhs_sample_idx], population_genotype_posteriors[rhs_sample_idx]);
                kl_divergences[lhs_sample_idx][rhs_sample_idx] = kl;
                kl_divergences[rhs_sample_idx][lhs_sample_idx] = kl;
            }
        }
        
        std::vector<std::vector<SampleName>> groups(num_groups);
        for (auto& group : groups) group.reserve(samples_.size());
        
        // Pick cluster centres
        std::size_t max_kl_sample1 {}, max_kl_sample2 {};
        double max_kl {0};
        for (std::size_t i {0}; i < num_groups; ++i) {
            for (std::size_t j {0}; j < i; ++j) {
                if (kl_divergences[i][j] > max_kl) {
                    max_kl_sample1 = i;
                    max_kl_sample2 = j;
                    max_kl = kl_divergences[i][j];
                }
            }
        }
        
        std::vector<bool> unassigned_samples(samples_.size(), false);
        
        groups[0].push_back(samples_[max_kl_sample1]);
        groups[1].push_back(samples_[max_kl_sample2]);
        
        unassigned_samples[max_kl_sample1] = true;
        unassigned_samples[max_kl_sample2] = true;
        
        for (std::size_t n {2}; n <= num_groups; ++n) {
            // TODO
        }
        
        for (std::size_t sample_idx {0}; sample_idx < samples_.size(); ++sample_idx) {
            if (!unassigned_samples[sample_idx]) {
                std::size_t best_group_idx {};
                double best_group_kl {max_kl};
                for (std::size_t j {0}; j < samples_.size(); ++j) {
                    if (j != sample_idx && unassigned_samples[j] && kl_divergences[sample_idx][j] < best_group_kl) {
                        best_group_idx = j;
                        best_group_kl = kl_divergences[sample_idx][j];
                    }
                }
                for (auto& group : groups) {
                    if (std::find(std::cbegin(group), std::cend(group), samples_[best_group_idx]) != std::cend(group)) {
                        group.push_back(samples_[sample_idx]);
                        break;
                    }
                }
            }
        }
        
        IndividualModel individual_model {prior_model_.germline_prior_model()};
        std::vector<ProbabilityVector> cluster_marginal_genotype_posteriors {};
        cluster_marginal_genotype_posteriors.reserve(num_groups);
        const auto haplotypes = extract_unique_elements(genotypes);
        for (const auto& group : groups) {
            const auto pooled_likelihoods = pool_likelihood(group, haplotypes, haplotype_likelihoods);
            auto cluster_inferences = individual_model.evaluate(genotypes, pooled_likelihoods);
            cluster_marginal_genotype_posteriors.push_back(std::move(cluster_inferences.posteriors.genotype_probabilities));
        }
        
        auto result = select_top_k_tuples(cluster_marginal_genotype_posteriors, config_.max_genotype_combinations);
        
        std::vector<int> counts(genotypes.size());
        result.erase(std::remove_if(std::begin(result), std::end(result), [&counts] (const auto& indices) {
            std::fill(std::begin(counts), std::end(counts), 0);
            for (auto idx : indices) ++counts[idx];
            return std::any_of(std::cbegin(counts), std::cend(counts), [] (auto count) { return count > 1; });
        }), std::end(result));
        
        return result;
    }
}

SingleCellModel::GenotypeCombinationVector
SingleCellModel::propose_all_genotype_combinations(const std::vector<Genotype<Haplotype>>& genotypes) const
{
    const auto num_groups = prior_model_.phylogeny().size();
    GenotypeCombinationVector result {};
    result.reserve(num_combinations(genotypes.size(), num_groups));
    GenotypeCombination tmp(num_groups);
    std::vector<bool> v(genotypes.size() * num_groups), hits(genotypes.size(), false);
    std::fill(std::begin(v), std::next(std::begin(v), num_groups), true);
    do {
        bool good {true};
        for (std::size_t i {0}, k {0}; k < num_groups; ++i) {
            if (v[i]) {
                if (i / genotypes.size() == k) {
                    tmp[k++] = i - genotypes.size() * (i / genotypes.size());
                } else {
                    good = false;
                    k = num_groups;
                }
            }
        }
        if (good) {
            for (auto idx : tmp) {
                if (hits[idx]) {
                    good = false;
                    break;
                } else {
                    hits[idx] = true;
                }
            }
            std::fill(std::begin(hits), std::end(hits), false);
            if (good) result.push_back(tmp);
        }
    } while (std::prev_permutation(std::begin(v), std::end(v)));
    return result;
}

VariationalBayesMixtureMixtureModel::LogProbabilityVector
SingleCellModel::calculate_genotype_priors(const GenotypeCombinationVector& genotype_combinations,
                                           const std::vector<Genotype<Haplotype>>& genotypes) const
{
    LogProbabilityVector result(genotype_combinations.size());
    std::transform(std::cbegin(genotype_combinations), std::cend(genotype_combinations), std::begin(result),
                   [&] (const auto& combination) {
                       using GenotypeRef = SingleCellPriorModel::GenotypeReference;
                       std::vector<GenotypeRef> genotype_refs {};
                       genotype_refs.reserve(combination.size());
                       std::transform(std::cbegin(combination), std::cend(combination), std::back_inserter(genotype_refs),
                                      [&] (auto idx) -> GenotypeRef { return genotypes[idx]; });
                       return prior_model_.evaluate(genotype_refs);
    });
    return result;
}

SingleCellModel::VBLikelihoodMatrix
SingleCellModel::make_likelihood_matrix(const GenotypeCombinationVector& genotype_combinations,
                                        const std::vector<Genotype<Haplotype>>& genotypes,
                                        const HaplotypeLikelihoodArray& haplotype_likelihoods) const
{
    VBLikelihoodMatrix result {};
    result.reserve(samples_.size());
    for (const auto& sample : samples_) {
        haplotype_likelihoods.prime(sample);
        VariationalBayesMixtureMixtureModel::GenotypeCombinationLikelihoodVector vb_combination_likelihoods {};
        vb_combination_likelihoods.reserve(genotype_combinations.size());
        for (const auto& genotype_combination : genotype_combinations) {
            VariationalBayesMixtureMixtureModel::GenotypeLikelihoodVector vb_genotype_likelihoods {};
            vb_genotype_likelihoods.reserve(genotype_combination.size());
            for (const auto genotype_idx : genotype_combination) {
                VariationalBayesMixtureMixtureModel::HaplotypeLikelihoodVector vb_haplotype_likelihoods {};
                vb_haplotype_likelihoods.reserve(genotypes[genotype_idx].ploidy());
                for (const auto& haplotype : genotypes[genotype_idx]) {
                    vb_haplotype_likelihoods.emplace_back(haplotype_likelihoods[haplotype]);
                }
                vb_genotype_likelihoods.push_back(std::move(vb_haplotype_likelihoods));
            }
            vb_combination_likelihoods.push_back(std::move(vb_genotype_likelihoods));
        }
        result.push_back(std::move(vb_combination_likelihoods));
    }
    return result;
}

namespace {

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

} // namespace

SingleCellModel::VBSeedVector
SingleCellModel::propose_seeds(const GenotypeCombinationVector& genotype_combinations,
                               const std::vector<Genotype<Haplotype>>& genotypes,
                               const VariationalBayesMixtureMixtureModel::LogProbabilityVector& genotype_combination_priors,
                               const HaplotypeLikelihoodArray& haplotype_likelihoods) const
{
    VBSeedVector result {};
    result.push_back(log_uniform_dist(genotype_combinations.size()));
    const auto k = std::min(std::size_t {config_.max_seeds}, genotype_combinations.size());
    std::vector<std::size_t> top_indices(k);
    std::iota(std::begin(top_indices), std::end(top_indices), 0u);
    make_point_seeds(genotype_combinations.size(), top_indices, result);
    return result;
}

} // namespace model
} // namespace octopus
