// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "single_cell_model.hpp"

#include <utility>
#include <cassert>

#include "subclone_model.hpp"

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

SingleCellModel::GenotypeCombinationVector
SingleCellModel::propose_genotype_combinations(const std::vector<Genotype<Haplotype>>& genotypes,
                                               const HaplotypeLikelihoodArray& haplotype_likelihoods) const
{
    const auto combination_size = prior_model_.phylogeny().size();
    GenotypeCombinationVector result {};
    result.reserve(config_.max_genotype_combinations);
    
    GenotypeCombination tmp(combination_size);
    std::vector<bool> v(genotypes.size() * combination_size), hits(genotypes.size(), false);
    std::fill(std::begin(v), std::next(std::begin(v), combination_size), true);
    do {
        bool good {true};
        for (std::size_t i {0}, k {0}; k < combination_size; ++i) {
            if (v[i]) {
                if (i / genotypes.size() == k) {
                    tmp[k++] = i - genotypes.size() * (i / genotypes.size());
                } else {
                    good = false;
                    k = combination_size;
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

SingleCellModel::VBSeedVector
SingleCellModel::propose_seeds(const GenotypeCombinationVector& genotype_combinations,
                               const std::vector<Genotype<Haplotype>>& genotypes,
                               const VariationalBayesMixtureMixtureModel::LogProbabilityVector& genotype_combination_priors,
                               const HaplotypeLikelihoodArray& haplotype_likelihoods) const
{
    VBSeedVector result {};
    result.emplace_back(genotype_combinations.size(), -std::log(genotype_combinations.size())); // Uniform
    
    // 1 5
//    for (std::size_t g_idx {0}; g_idx < genotype_combinations.size(); ++g_idx) {
//        if ((genotype_combinations[g_idx][0] == 1 && genotype_combinations[g_idx][1] == 5)) {
//            result.emplace_back(genotype_combinations.size(), std::log((1.0 - 0.999999) / (genotype_combinations.size() - 1)));
//            result.back()[g_idx] = 0.999999;
//        }
//    }
    
    return result;
}

} // namespace model
} // namespace octopus
