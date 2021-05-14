// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "family_model.hpp"

#include "utils/maths.hpp"
#include "utils/select_top_k.hpp"
#include "utils/append.hpp"
#include "constant_mixture_genotype_likelihood_model.hpp"
#include "population_model.hpp"

namespace octopus { namespace model {

FamilyModel::FamilyModel(const Pedigree& family,
                         const PedigreePriorModel& prior_model,
                         Options options,
                         boost::optional<logging::DebugLogger> debug_log)
: family_ {family}
, prior_model_ {prior_model}
, options_ {options}
, debug_log_ {debug_log}
{}

const PedigreePriorModel& FamilyModel::prior_model() const noexcept
{
    return prior_model_;
}

namespace {

using LogProbability = double;
using GenotypeLogLikelihoodVector  = std::vector<LogProbability>;
using GenotypeLogLikelihoodMatrix  = std::vector<GenotypeLogLikelihoodVector>;

GenotypeLogLikelihoodMatrix
compute_genotype_log_likelihoods(const std::vector<SampleName>& samples,
                                 const FamilyModel::GenotypeVector& genotypes,
                                 const HaplotypeLikelihoodArray& haplotype_likelihoods)
{
    assert(!genotypes.empty());
    ConstantMixtureGenotypeLikelihoodModel likelihood_model {haplotype_likelihoods};
    GenotypeLogLikelihoodMatrix result {};
    result.reserve(samples.size());
    std::transform(std::cbegin(samples), std::cend(samples), std::back_inserter(result), [&] (const auto& sample) {
        GenotypeLogLikelihoodVector likelihoods(genotypes.size());
        haplotype_likelihoods.prime(sample);
        std::transform(std::cbegin(genotypes), std::cend(genotypes), std::begin(likelihoods),
                       [&] (const auto& genotype) { return likelihood_model.evaluate(genotype); });
        return likelihoods;
    });
    return result;
}

boost::optional<std::size_t> compute_num_combinations(const std::size_t num_genotypes, const std::size_t num_samples)
{
    if (maths::is_safe_ipow(num_genotypes, num_samples)) {
        return maths::ipow(num_genotypes, num_samples);
    } else {
        return boost::none;
    }
}

using GenotypeCombinationVector = std::vector<unsigned>;
using GenotypeCombinationMatrix = std::vector<GenotypeCombinationVector>;

auto generate_all_genotype_combinations(const std::size_t num_genotypes, const std::size_t num_samples)
{
    GenotypeCombinationMatrix result {};
    const auto num_combinations = compute_num_combinations(num_genotypes, num_samples);
    if (!num_combinations) throw std::overflow_error {"generate_all_genotype_combinations overflowed"};
    result.reserve(*num_combinations);
    GenotypeCombinationVector tmp(num_samples);
    std::vector<bool> v(num_genotypes * num_samples);
    std::fill(std::begin(v), std::next(std::begin(v), num_samples), true);
    do {
        bool good {true};
        for (std::size_t i {0}, k {0}; k < num_samples; ++i) {
            if (v[i]) {
                if (i / num_genotypes == k) {
                    tmp[k++] = i - num_genotypes * (i / num_genotypes);
                } else {
                    good = false;
                    k = num_samples;
                }
            }
        }
        if (good) result.push_back(tmp);
    } while (std::prev_permutation(std::begin(v), std::end(v)));
    return result;
}

template <typename Range>
boost::optional<std::size_t> find_hom_ref_idx(const Range& genotypes)
{
    auto itr = std::find_if(std::cbegin(genotypes), std::cend(genotypes),
                            [] (const auto& g) { return is_homozygous_reference(g); });
    if (itr != std::cend(genotypes)) {
        return std::distance(std::cbegin(genotypes), itr);
    } else {
        return boost::none;
    }
}

template <typename Container>
auto sum(const Container& values)
{
    return std::accumulate(std::cbegin(values), std::cend(values), 0.0);
}

void fill(const GenotypeLogLikelihoodMatrix& genotype_likelihoods,
          const GenotypeCombinationVector& combination,
          GenotypeLogLikelihoodVector& result)
{
    assert(result.size() == combination.size());
    for (std::size_t s {0}; s < combination.size(); ++s) {
        result[s] = genotype_likelihoods[s][combination[s]];
    }
}

} // namespace

FamilyModel::InferredLatents
FamilyModel::evaluate(const SampleVector& samples,
                      const MappableBlock<Haplotype>& haplotypes,
                      const GenotypeVector& genotypes,
                      const HaplotypeLikelihoodArray& haplotype_likelihoods) const
{
    assert(!genotypes.empty());
    const auto max_possible_genotype_combinations = compute_num_combinations(genotypes.size(), samples.size());
    GenotypeLogLikelihoodMatrix genotype_likelihoods {};
    GenotypeCombinationMatrix genotype_combinations {};
    if (!max_possible_genotype_combinations || *max_possible_genotype_combinations < options_.max_genotype_combinations) {
        genotype_likelihoods = compute_genotype_log_likelihoods(samples, genotypes, haplotype_likelihoods);
        genotype_combinations = generate_all_genotype_combinations(genotypes.size(), samples.size());
    } else {
        PopulationModel::Options population_model_options {};
        population_model_options.max_genotype_combinations = options_.max_genotype_combinations;
        PopulationModel population_model {prior_model_.population_model(), population_model_options, debug_log_};
        auto population_inferences = population_model.evaluate(samples, haplotypes, genotypes, haplotype_likelihoods);
        genotype_combinations = select_top_k_tuples<unsigned>(population_inferences.posteriors.marginal_genotype_probabilities, *options_.max_genotype_combinations);
        const auto hom_ref_genotype_idx = find_hom_ref_idx(genotypes);
        if (hom_ref_genotype_idx) {
            // Add some reference combinations to get nicer QUALs
            std::vector<GenotypeCombinationVector> reference_combinations {};
            GenotypeCombinationVector ref_indices(samples.size(), *hom_ref_genotype_idx);
            if (std::find(std::cbegin(genotype_combinations), std::cend(genotype_combinations), ref_indices) == std::cend(genotype_combinations)) {
                reference_combinations.push_back(std::move(ref_indices));
            }
            for (std::size_t s {0}; s < samples.size(); ++s) {
                const auto is_sample_homref = [=] (const auto& combination) { return combination[s] == *hom_ref_genotype_idx; };
                if (std::find_if(std::cbegin(genotype_combinations), std::cend(genotype_combinations), is_sample_homref) == std::cend(genotype_combinations)) {
                    auto new_combination = genotype_combinations.front(); // the best
                    new_combination[s] = *hom_ref_genotype_idx;
                    if (std::find(std::cbegin(reference_combinations), std::cend(reference_combinations), new_combination) == std::cend(reference_combinations)) {
                        reference_combinations.push_back(std::move(new_combination));
                    }
                }
            }
            utils::append(std::move(reference_combinations), genotype_combinations);
        }
        genotype_likelihoods = std::move(population_inferences.posteriors.genotype_log_likelihoods);
    }
    std::vector<double> posteriors(genotype_combinations.size());
    GenotypeLogLikelihoodVector likelihoods_buffer(genotype_likelihoods.size());
    for (std::size_t g {0}; g < genotype_combinations.size(); ++g) {
        fill(genotype_likelihoods, genotype_combinations[g], likelihoods_buffer);
        std::vector<Genotype<IndexedHaplotype<>>> combination {};
        combination.reserve(samples.size());
        for (auto i : genotype_combinations[g]) {
            combination.push_back(genotypes[i]);
        }
        posteriors[g] = prior_model_.evaluate(combination) + sum(likelihoods_buffer);
    }
    InferredLatents result {};
    result.log_evidence = maths::normalise_logs(posteriors);
    result.posteriors.joint_genotype_probabilities.reserve(genotype_combinations.size());
    for (std::size_t g {0}; g < genotype_combinations.size(); ++g) {
        result.posteriors.joint_genotype_probabilities.push_back({genotype_combinations[g], posteriors[g]});
    }
    return result;
}

FamilyModel::InferredLatents
FamilyModel::evaluate(const SampleVector& samples,
                      const std::vector<unsigned>& sample_ploidies,
                      const MappableBlock<Haplotype>& haplotypes,
                      const GenotypeVector& genotypes,
                      const HaplotypeLikelihoodArray& haplotype_likelihoods) const
{
    return {};
}

} // namespace model
} // namespace octopus
