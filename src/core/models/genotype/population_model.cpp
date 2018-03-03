// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "population_model.hpp"

#include <unordered_map>
#include <utility>
#include <algorithm>
#include <cmath>
#include <limits>
#include <cassert>

#include "utils/maths.hpp"
#include "germline_likelihood_model.hpp"

namespace octopus { namespace model {

PopulationModel::PopulationModel(const PopulationPriorModel& prior_model,
                                 boost::optional<logging::DebugLogger> debug_log)
: prior_model_ {prior_model}
, debug_log_ {debug_log}
{}

PopulationModel::PopulationModel(const PopulationPriorModel& prior_model,
                                 Options options,
                                 boost::optional<logging::DebugLogger> debug_log)
: options_ {options}
, prior_model_ {prior_model}
, debug_log_ {debug_log}
{}

const PopulationPriorModel& PopulationModel::prior_model() const noexcept
{
    return prior_model_;
}

namespace {

using GenotypeLogLikelihoodVector  = std::vector<double>;
using GenotypeLogLikelihoodMatrix  = std::vector<GenotypeLogLikelihoodVector>;

struct GenotypeLogProbability
{
    const Genotype<Haplotype>& genotype;
    double log_probability;
};
using GenotypeLogMarginalVector = std::vector<GenotypeLogProbability>;

using GenotypeMarginalPosteriorVector  = std::vector<double>;
using GenotypeMarginalPosteriorMatrix  = std::vector<GenotypeMarginalPosteriorVector>; // for each sample

using InverseGenotypeTable = std::vector<std::vector<std::size_t>>;

auto make_inverse_genotype_table(const std::vector<Haplotype>& haplotypes,
                                 const std::vector<Genotype<Haplotype>>& genotypes)
{
    assert(!haplotypes.empty() && !genotypes.empty());
    using HaplotypeReference = std::reference_wrapper<const Haplotype>;
    std::unordered_map<HaplotypeReference, std::vector<std::size_t>> result_map {haplotypes.size()};
    const auto cardinality = element_cardinality_in_genotypes(static_cast<unsigned>(haplotypes.size()),
                                                              genotypes.front().ploidy());
    for (const auto& haplotype : haplotypes) {
        result_map[std::cref(haplotype)].reserve(cardinality);
    }
    for (std::size_t i {0}; i < genotypes.size(); ++i) {
        for (const auto& haplotype : genotypes[i]) {
            result_map.at(haplotype).emplace_back(i);
        }
    }
    InverseGenotypeTable result {};
    result.reserve(haplotypes.size());
    for (const auto& haplotype : haplotypes) {
        result.emplace_back(std::move(result_map.at(haplotype)));
    }
    return result;
}

auto make_inverse_genotype_table(const std::vector<std::vector<unsigned>>& genotype_indices,
                                 const std::size_t num_haplotypes)
{
    InverseGenotypeTable result(num_haplotypes);
    for (const auto& genotype : genotype_indices) {
        for (auto idx : genotype) {
            if (result[idx].empty() || result[idx].back() != idx) {
                result[idx].push_back(idx);
            }
        }
    }
    for (auto& entry : result) entry.shrink_to_fit();
    return result;
}

using HaplotypeReference    = std::reference_wrapper<const Haplotype>;
using HaplotypeFrequencyMap = std::unordered_map<HaplotypeReference, double>;

double calculate_frequency_update_norm(const std::size_t num_samples, const unsigned ploidy)
{
    return static_cast<double>(num_samples) * ploidy;
}

struct EMOptions
{
    unsigned max_iterations;
    double epsilon;
};

struct ModelConstants
{
    const std::vector<Haplotype>& haplotypes;
    const std::vector<Genotype<Haplotype>>& genotypes;
    const GenotypeLogLikelihoodMatrix& genotype_log_likilhoods;
    const unsigned ploidy;
    const double frequency_update_norm;
    const InverseGenotypeTable genotypes_containing_haplotypes;
    
    ModelConstants(const std::vector<Haplotype>& haplotypes,
                   const std::vector<Genotype<Haplotype>>& genotypes,
                   const GenotypeLogLikelihoodMatrix& genotype_log_likilhoods)
    : haplotypes {haplotypes}
    , genotypes {genotypes}
    , genotype_log_likilhoods {genotype_log_likilhoods}
    , ploidy {genotypes.front().ploidy()}
    , frequency_update_norm {calculate_frequency_update_norm(genotype_log_likilhoods.size(), ploidy)}
    , genotypes_containing_haplotypes {make_inverse_genotype_table(haplotypes, genotypes)}
    {}
    ModelConstants(const std::vector<Haplotype>& haplotypes,
                   const std::vector<Genotype<Haplotype>>& genotypes,
                   const std::vector<std::vector<unsigned>>& genotype_indices,
                   const GenotypeLogLikelihoodMatrix& genotype_log_likilhoods)
    : haplotypes {haplotypes}
    , genotypes {genotypes}
    , genotype_log_likilhoods {genotype_log_likilhoods}
    , ploidy {genotypes.front().ploidy()}
    , frequency_update_norm {calculate_frequency_update_norm(genotype_log_likilhoods.size(), ploidy)}
    , genotypes_containing_haplotypes {make_inverse_genotype_table(genotype_indices, haplotypes.size())}
    {}
};

HaplotypeFrequencyMap
init_haplotype_frequencies(const ModelConstants& constants)
{
    HaplotypeFrequencyMap result {constants.haplotypes.size()};
    for (const auto& haplotype : constants.haplotypes) {
        result.emplace(haplotype, 1.0 / constants.haplotypes.size());
    }
    return result;
}

template <typename Genotype, typename Map>
double log_hardy_weinberg_haploid(const Genotype& genotype,
                                  const Map& haplotype_frequencies)
{
    return std::log(haplotype_frequencies.at(genotype[0]));
}

template <typename Genotype, typename Map>
double log_hardy_weinberg_diploid(const Genotype& genotype,
                                  const Map& haplotype_frequencies)
{
    if (genotype.is_homozygous()) {
        return 2 * std::log(haplotype_frequencies.at(genotype[0]));
    }
    static const double ln2 {std::log(2.0)};
    return std::log(haplotype_frequencies.at(genotype[0])) + std::log(haplotype_frequencies.at(genotype[1])) + ln2;
}

template <typename Genotype, typename Map>
double log_hardy_weinberg_triploid(const Genotype& genotype,
                                   const Map& haplotype_frequencies)
{
    // TODO: optimise this case
    auto unique_haplotypes = genotype.copy_unique();
    std::vector<unsigned> occurences {};
    occurences.reserve(unique_haplotypes.size());
    double r {0};
    for (const auto& haplotype : unique_haplotypes) {
        auto num_occurences = genotype.count(haplotype);
        occurences.push_back(num_occurences);
        r += num_occurences * std::log(haplotype_frequencies.at(haplotype));
    }
    return maths::log_multinomial_coefficient<double>(occurences) + r;
}

template <typename Genotype, typename Map>
double log_hardy_weinberg_polyploid(const Genotype& genotype,
                                    const Map& haplotype_frequencies)
{
    auto unique_haplotypes = genotype.copy_unique();
    std::vector<unsigned> occurences {};
    occurences.reserve(unique_haplotypes.size());
    double r {0};
    for (const auto& haplotype : unique_haplotypes) {
        auto num_occurences = genotype.count(haplotype);
        occurences.push_back(num_occurences);
        r += num_occurences * std::log(haplotype_frequencies.at(haplotype));
    }
    return maths::log_multinomial_coefficient<double>(occurences) + r;
}

// TODO: improve this, possible bottleneck in EM update at the moment
template <typename Genotype, typename Map>
double log_hardy_weinberg(const Genotype& genotype, const Map& haplotype_frequencies)
{
    switch (genotype.ploidy()) {
        case 1 : return log_hardy_weinberg_haploid(genotype, haplotype_frequencies);
        case 2 : return log_hardy_weinberg_diploid(genotype, haplotype_frequencies);
        case 3 : return log_hardy_weinberg_triploid(genotype, haplotype_frequencies);
        default: return log_hardy_weinberg_polyploid(genotype, haplotype_frequencies);
    }
}

GenotypeLogLikelihoodMatrix
compute_genotype_log_likelihoods(const std::vector<SampleName>& samples,
                                 const std::vector<Genotype<Haplotype>>& genotypes,
                                 const HaplotypeLikelihoodCache& haplotype_likelihoods)
{
    assert(!genotypes.empty());
    GermlineLikelihoodModel likelihood_model {haplotype_likelihoods};
    GenotypeLogLikelihoodMatrix result {};
    result.reserve(samples.size());
    std::transform(std::cbegin(samples), std::cend(samples), std::back_inserter(result),
                   [&genotypes, &haplotype_likelihoods, &likelihood_model] (const auto& sample) {
                       GenotypeLogLikelihoodVector likelihoods(genotypes.size());
                       haplotype_likelihoods.prime(sample);
                       std::transform(std::cbegin(genotypes), std::cend(genotypes), std::begin(likelihoods),
                                      [&likelihood_model] (const auto& genotype) {
                                          return likelihood_model.evaluate(genotype);
                                      });
                       return likelihoods;
                   });
    return result;
}

GenotypeLogMarginalVector
init_genotype_log_marginals(const std::vector<Genotype<Haplotype>>& genotypes,
                            const HaplotypeFrequencyMap& haplotype_frequencies)
{
    GenotypeLogMarginalVector result {};
    result.reserve(genotypes.size());
    for (const auto& genotype : genotypes) {
        result.push_back({genotype, log_hardy_weinberg(genotype, haplotype_frequencies)});
    }
    return result;
}

void update_genotype_log_marginals(GenotypeLogMarginalVector& current_log_marginals,
                                   const HaplotypeFrequencyMap& haplotype_frequencies)
{
    std::for_each(std::begin(current_log_marginals), std::end(current_log_marginals),
                  [&haplotype_frequencies] (auto& p) {
                      p.log_probability = log_hardy_weinberg(p.genotype, haplotype_frequencies);
                  });
}

GenotypeMarginalPosteriorMatrix
init_genotype_posteriors(const GenotypeLogMarginalVector& genotype_log_marginals,
                         const GenotypeLogLikelihoodMatrix& genotype_log_likilhoods)
{
    GenotypeMarginalPosteriorMatrix result {};
    result.reserve(genotype_log_likilhoods.size());
    for (const auto& sample_genotype_log_likilhoods : genotype_log_likilhoods) {
        GenotypeMarginalPosteriorVector posteriors(genotype_log_marginals.size());
        std::transform(std::cbegin(genotype_log_marginals), std::cend(genotype_log_marginals),
                       std::cbegin(sample_genotype_log_likilhoods),
                       std::begin(posteriors),
                       [] (const auto& genotype_log_marginal, const auto genotype_log_likilhood) {
                           return genotype_log_marginal.log_probability + genotype_log_likilhood;
                       });
        maths::normalise_exp(posteriors);
        result.emplace_back(std::move(posteriors));
    }
    return result;
}

void update_genotype_posteriors(GenotypeMarginalPosteriorMatrix& current_genotype_posteriors,
                                const GenotypeLogMarginalVector& genotype_log_marginals,
                                const GenotypeLogLikelihoodMatrix& genotype_log_likilhoods)
{
    auto it = std::cbegin(genotype_log_likilhoods);
    for (auto& sample_genotype_posteriors : current_genotype_posteriors) {
        std::transform(std::cbegin(genotype_log_marginals), std::cend(genotype_log_marginals),
                       std::cbegin(*it++),
                       std::begin(sample_genotype_posteriors),
                       [] (const auto& log_marginal, const auto& log_likilhood) {
                           return log_marginal.log_probability + log_likilhood;
                       });
        maths::normalise_exp(sample_genotype_posteriors);
    }
}

auto collapse_genotype_posteriors(const GenotypeMarginalPosteriorMatrix& genotype_posteriors)
{
    assert(!genotype_posteriors.empty());
    std::vector<double> result(genotype_posteriors.front().size());
    for (const auto& sample_posteriors : genotype_posteriors) {
        std::transform(std::cbegin(result), std::cend(result), std::cbegin(sample_posteriors), std::begin(result),
                       [] (const auto curr, const auto p) { return curr + p; });
    }
    return result;
}

double update_haplotype_frequencies(const std::vector<Haplotype>& haplotypes,
                                    HaplotypeFrequencyMap& current_haplotype_frequencies,
                                    const GenotypeMarginalPosteriorMatrix& genotype_posteriors,
                                    const InverseGenotypeTable& genotypes_containing_haplotypes,
                                    const double frequency_update_norm)
{
    const auto collaped_posteriors = collapse_genotype_posteriors(genotype_posteriors);
    double max_frequency_change {0};
    for (std::size_t i {0}; i < haplotypes.size(); ++i) {
        auto& current_frequency = current_haplotype_frequencies.at(haplotypes[i]);
        double new_frequency {0};
        for (const auto& genotype_index : genotypes_containing_haplotypes[i]) {
            new_frequency += collaped_posteriors[genotype_index];
        }
        new_frequency /= frequency_update_norm;
        const auto frequency_change = std::abs(current_frequency - new_frequency);
        if (frequency_change > max_frequency_change) {
            max_frequency_change = frequency_change;
        }
        current_frequency = new_frequency;
    }
    return max_frequency_change;
}

double do_em_iteration(GenotypeMarginalPosteriorMatrix& genotype_posteriors,
                       HaplotypeFrequencyMap& haplotype_frequencies,
                       GenotypeLogMarginalVector& genotype_log_marginals,
                       const ModelConstants& constants)
{
    const auto max_change = update_haplotype_frequencies(constants.haplotypes,
                                                         haplotype_frequencies,
                                                         genotype_posteriors,
                                                         constants.genotypes_containing_haplotypes,
                                                         constants.frequency_update_norm);
    update_genotype_log_marginals(genotype_log_marginals, haplotype_frequencies);
    update_genotype_posteriors(genotype_posteriors, genotype_log_marginals,
                               constants.genotype_log_likilhoods);
    return max_change;
}

void run_em(GenotypeMarginalPosteriorMatrix& genotype_posteriors,
            HaplotypeFrequencyMap& haplotype_frequencies,
            GenotypeLogMarginalVector& genotype_log_marginals,
            const ModelConstants& constants, const EMOptions options,
            boost::optional<logging::TraceLogger> trace_log = boost::none)
{
    for (unsigned n {1}; n <= options.max_iterations; ++n) {
        const auto max_change = do_em_iteration(genotype_posteriors, haplotype_frequencies,
                                                genotype_log_marginals,constants);
        if (max_change <= options.epsilon) break;
    }
}

auto compute_approx_genotype_marginal_posteriors(const std::vector<Haplotype>& haplotypes,
                                                 const std::vector<Genotype<Haplotype>>& genotypes,
                                                 const GenotypeLogLikelihoodMatrix& genotype_likelihoods,
                                                 const EMOptions options)
{
    const ModelConstants constants {haplotypes, genotypes, genotype_likelihoods};
    auto haplotype_frequencies = init_haplotype_frequencies(constants);
    auto genotype_log_marginals = init_genotype_log_marginals(genotypes, haplotype_frequencies);
    auto result = init_genotype_posteriors(genotype_log_marginals, genotype_likelihoods);
    run_em(result, haplotype_frequencies, genotype_log_marginals, constants, options);
    return result;
}

auto compute_approx_genotype_marginal_posteriors(const std::vector<Genotype<Haplotype>>& genotypes,
                                                 const GenotypeLogLikelihoodMatrix& genotype_likelihoods,
                                                 const EMOptions options)
{
    const auto haplotypes = extract_unique_elements(genotypes);
    return compute_approx_genotype_marginal_posteriors(haplotypes, genotypes, genotype_likelihoods, options);
}

using GenotypeCombinationVector = std::vector<std::size_t>;
using GenotypeCombinationMatrix = std::vector<GenotypeCombinationVector>;

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

auto generate_all_genotype_combinations(const std::size_t num_genotypes, const std::size_t num_samples)
{
    GenotypeCombinationMatrix result {};
    result.reserve(num_combinations(num_genotypes, num_samples));
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

template <typename T>
auto index(const std::vector<T>& values)
{
    std::vector<std::pair<T, std::size_t>> result(values.size());
    for (std::size_t i {0}; i < values.size(); ++i) {
        result[i] = std::make_pair(values[i], i);
    }
    return result;
}

auto index_and_sort(const GenotypeMarginalPosteriorVector& genotype_posteriors)
{
    auto result = index(genotype_posteriors);
    std::sort(std::begin(result), std::end(result));
    return result;
}

auto index_and_sort(const GenotypeMarginalPosteriorMatrix& marginals)
{
    std::vector<std::vector<std::pair<double, std::size_t>>> result {};
    result.reserve(marginals.size());
    for (const auto& sample_marginals : marginals) {
        result.push_back(index_and_sort(sample_marginals));
    }
    return result;
}

auto propose_joint_genotypes(const std::size_t num_genotypes,
                             const GenotypeMarginalPosteriorMatrix& em_genotype_marginals,
                             const std::size_t max_joint_genotypes)
{
    const auto num_samples = em_genotype_marginals.size();
    assert(max_joint_genotypes >= num_samples * num_genotypes);
    const auto num_joint_genotypes = num_combinations(num_genotypes, num_samples);
    if (num_joint_genotypes <= max_joint_genotypes) {
        return generate_all_genotype_combinations(num_genotypes, num_samples);
    }
    const auto ranked_marginals = index_and_sort(em_genotype_marginals);
    GenotypeCombinationMatrix result {};
    result.reserve(max_joint_genotypes);
    auto remaining_combinations = max_joint_genotypes;
    std::vector<std::size_t> map_indepedent_joint(num_samples);
    for (std::size_t s {0}; s < num_samples; ++s) {
        map_indepedent_joint[s] = ranked_marginals[s].front().second;
    }
    // Ensure each genotype is represented at least once for each sample
    for (std::size_t s {0}; s < num_samples; ++s) {
        for (std::size_t g = s > 0 ? 1 : 0; g < num_genotypes; ++g) {
            result.push_back(map_indepedent_joint);
            result.back()[s] = g;
            --remaining_combinations;
        }
    }
    if (remaining_combinations > 0) {
//        std::generate_n(std::back_inserter(result), remaining_combinations, [&] () {
//
//        });
    }
    return result;
}

template <typename Container>
auto sum(const Container& values)
{
    return std::accumulate(std::cbegin(values), std::cend(values), 0.0);
}

void fill(const GenotypeLogLikelihoodMatrix& genotype_likelihoods,
          const GenotypeCombinationVector& indices,
          GenotypeLogLikelihoodVector& result)
{
    assert(result.size() == indices.size());
    for (std::size_t s {0}; s < indices.size(); ++s) {
        result[s] = genotype_likelihoods[s][indices[s]];
    }
}

using GenotypeReferenceVector = std::vector<std::reference_wrapper<const Genotype<Haplotype>>>;

void fill(const std::vector<Genotype<Haplotype>>& genotypes,
          const GenotypeCombinationVector& indices,
          GenotypeReferenceVector& result)
{
    result.clear();
    std::transform(std::cbegin(indices), std::cend(indices), std::back_inserter(result),
                   [&genotypes] (const auto index) { return std::cref(genotypes[index]); });
}

auto calculate_posteriors(const std::vector<Genotype<Haplotype>>& genotypes,
                          const GenotypeCombinationMatrix& joint_genotypes,
                          const GenotypeLogLikelihoodMatrix& genotype_likelihoods,
                          const PopulationPriorModel& prior_model)
{
    std::vector<double> result {};
    GenotypeLogLikelihoodVector likelihoods_buffer(genotype_likelihoods.size());
    GenotypeReferenceVector genotypes_refs {};
    for (const auto& indices : joint_genotypes) {
        fill(genotype_likelihoods, indices, likelihoods_buffer);
        fill(genotypes, indices, genotypes_refs);
        result.push_back(prior_model.evaluate(genotypes_refs) + sum(likelihoods_buffer));
    }
    const auto norm = maths::normalise_exp(result);
    return std::make_pair(std::move(result), norm);
}

void calculate_posterior_marginals(const std::vector<Genotype<Haplotype>>& genotypes,
                                   const GenotypeCombinationMatrix& joint_genotypes,
                                   const GenotypeLogLikelihoodMatrix& genotype_likelihoods,
                                   const PopulationPriorModel& prior_model,
                                   PopulationModel::InferredLatents& result)
{
    std::vector<double> joint_posteriors; double norm;
    std::tie(joint_posteriors, norm) = calculate_posteriors(genotypes, joint_genotypes, genotype_likelihoods, prior_model);
    assert(joint_posteriors.size() == joint_genotypes.size());
    const auto num_samples = genotype_likelihoods.size();
    std::vector<std::vector<double>> marginals(num_samples, std::vector<double>(genotypes.size(), 0.0));
    for (std::size_t i {0}; i < joint_genotypes.size(); ++i) {
        assert(joint_genotypes[i].size() == num_samples);
        for (std::size_t s {0}; s < num_samples; ++s) {
            marginals[s][joint_genotypes[i][s]] += joint_posteriors[i];
        }
    }
    result.posteriors.marginal_genotype_probabilities = std::move(marginals);
    result.log_evidence = norm;
}

} // namespace

PopulationModel::InferredLatents
PopulationModel::evaluate(const SampleVector& samples, const GenotypeVector& genotypes,
                          const HaplotypeLikelihoodCache& haplotype_likelihoods) const
{
    assert(!genotypes.empty());
    const auto genotype_log_likelihoods = compute_genotype_log_likelihoods(samples, genotypes, haplotype_likelihoods);
    const auto num_joint_genotypes = num_combinations(genotypes.size(), samples.size());
    InferredLatents result;
    if (num_joint_genotypes <= options_.max_joint_genotypes) {
        const auto joint_genotypes = generate_all_genotype_combinations(genotypes.size(), samples.size());
        calculate_posterior_marginals(genotypes, joint_genotypes, genotype_log_likelihoods, prior_model_, result);
    } else {
        const EMOptions em_options {options_.max_em_iterations, options_.em_epsilon};
        const auto em_genotype_marginals = compute_approx_genotype_marginal_posteriors(genotypes, genotype_log_likelihoods, em_options);
        const auto joint_genotypes = propose_joint_genotypes(genotypes.size(), em_genotype_marginals, options_.max_joint_genotypes);
        calculate_posterior_marginals(genotypes, joint_genotypes, genotype_log_likelihoods, prior_model_, result);
    }
    return result;
}

PopulationModel::InferredLatents
PopulationModel::evaluate(const SampleVector& samples,
                          const GenotypeVector& genotypes,
                          const std::vector<std::vector<unsigned>>& genotype_indices,
                          const HaplotypeLikelihoodCache& haplotype_likelihoods) const
{
    assert(!genotypes.empty());
    const auto genotype_log_likelihoods = compute_genotype_log_likelihoods(samples, genotypes, haplotype_likelihoods);
    InferredLatents result;
    return result;
}

PopulationModel::InferredLatents
PopulationModel::evaluate(const SampleVector& samples,
                          const std::vector<GenotypeVectorReference>& genotypes,
                          const HaplotypeLikelihoodCache& haplotype_likelihoods) const
{
    return InferredLatents {};
}

namespace debug {
    
} // namespace debug

} // namespace model
} // namespace octopus
