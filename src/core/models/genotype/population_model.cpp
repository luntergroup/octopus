// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "population_model.hpp"

#include <unordered_map>
#include <utility>
#include <algorithm>
#include <cmath>
#include <limits>
#include <cassert>
#include <exception>

#include "utils/maths.hpp"
#include "utils/select_top_k.hpp"
#include "utils/concat.hpp"
#include "constant_mixture_genotype_likelihood_model.hpp"
#include "hardy_weinberg_model.hpp"

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

using LogProbability = double;
using GenotypeLogLikelihoodVector  = std::vector<LogProbability>;
using GenotypeLogLikelihoodMatrix  = std::vector<GenotypeLogLikelihoodVector>;

struct GenotypeLogProbability
{
    const Genotype<IndexedHaplotype<>>& genotype;
    double log_probability;
};
using GenotypeLogMarginalVector = std::vector<GenotypeLogProbability>;

using GenotypeMarginalPosteriorVector  = std::vector<double>;
using GenotypeMarginalPosteriorMatrix  = std::vector<GenotypeMarginalPosteriorVector>; // for each sample

using InverseGenotypeTable = std::vector<std::vector<std::size_t>>;

auto make_inverse_genotype_table(const MappableBlock<Genotype<IndexedHaplotype<>>>& genotypes,
                                 const std::size_t num_haplotypes)
{
    InverseGenotypeTable result(num_haplotypes);
    for (auto& indices : result) indices.reserve(genotypes.size() / num_haplotypes);
    for (std::size_t genotype_idx {0}; genotype_idx < genotypes.size(); ++genotype_idx) {
        for (const auto& haplotype : genotypes[genotype_idx]) {
            result[index_of(haplotype)].push_back(genotype_idx);
        }
    }
    for (auto& indices : result) {
        std::sort(std::begin(indices), std::end(indices));
        indices.erase(std::unique(std::begin(indices), std::end(indices)), std::end(indices));
        indices.shrink_to_fit();
    }
    return result;
}

double calculate_frequency_update_norm(const std::size_t num_samples, const unsigned ploidy) noexcept
{
    return static_cast<double>(num_samples) * ploidy;
}

double calculate_frequency_update_norm(const std::vector<unsigned>& sample_ploidies) noexcept
{
    return std::accumulate(std::cbegin(sample_ploidies), std::cend(sample_ploidies), 1.0, std::multiplies<> {});
}

struct EMOptions
{
    unsigned max_iterations;
    double epsilon;
};

struct ModelConstants
{
    const PopulationModel::GenotypeVector& genotypes;
    const GenotypeLogLikelihoodMatrix& genotype_log_likilhoods;
    const InverseGenotypeTable genotypes_containing_haplotypes;
    const std::size_t num_haplotypes;
    const double frequency_update_norm;
    
    ModelConstants(const MappableBlock<Haplotype>& haplotypes,
                   const PopulationModel::GenotypeVector& genotypes,
                   const GenotypeLogLikelihoodMatrix& genotype_log_likilhoods)
    : genotypes {genotypes}
    , genotype_log_likilhoods {genotype_log_likilhoods}
    , genotypes_containing_haplotypes {make_inverse_genotype_table(genotypes, haplotypes.size())}
    , num_haplotypes {haplotypes.size()}
    , frequency_update_norm {calculate_frequency_update_norm(genotype_log_likilhoods.size(), genotypes.front().ploidy())}
    {}
    ModelConstants(const MappableBlock<Haplotype>& haplotypes,
                   const PopulationModel::GenotypeVector& genotypes,
                   const GenotypeLogLikelihoodMatrix& genotype_log_likilhoods,
                   const std::vector<unsigned>& sample_ploidies)
    : genotypes {genotypes}
    , genotype_log_likilhoods {genotype_log_likilhoods}
    , genotypes_containing_haplotypes {make_inverse_genotype_table(genotypes, haplotypes.size())}
    , num_haplotypes {haplotypes.size()}
    , frequency_update_norm {calculate_frequency_update_norm(sample_ploidies)}
    {}
};

HardyWeinbergModel make_hardy_weinberg_model(const ModelConstants& constants)
{
    HardyWeinbergModel::HaplotypeFrequencyVector frequencies(constants.num_haplotypes);
    for (std::size_t idx {0}; idx < constants.num_haplotypes; ++idx) {
        frequencies[idx] = 1.0 / constants.num_haplotypes;
    }
    return HardyWeinbergModel {std::move(frequencies)};
}

GenotypeLogLikelihoodMatrix
compute_genotype_log_likelihoods(const std::vector<SampleName>& samples,
                                 const PopulationModel::GenotypeVector& genotypes,
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

GenotypeLogLikelihoodMatrix
compute_genotype_log_likelihoods(const std::vector<SampleName>& samples,
                                 const PopulationModel::GenotypeVector& genotypes,
                                 const HaplotypeLikelihoodArray& haplotype_likelihoods,
                                 const std::vector<std::vector<bool>>& sample_genotype_masks)
{
    assert(!genotypes.empty());
    ConstantMixtureGenotypeLikelihoodModel likelihood_model {haplotype_likelihoods};
    GenotypeLogLikelihoodMatrix result {};
    result.reserve(samples.size());
    std::transform(std::cbegin(samples), std::cend(samples), std::cbegin(sample_genotype_masks),
                   std::back_inserter(result), [&] (const auto& sample, const auto& mask) {
        GenotypeLogLikelihoodVector likelihoods(genotypes.size());
        haplotype_likelihoods.prime(sample);
        std::transform(std::cbegin(genotypes), std::cend(genotypes), std::cbegin(mask), std::begin(likelihoods),
                       [&] (const auto& genotype, bool ok) {
                           return ok ? likelihood_model.evaluate(genotype) : std::numeric_limits<LogProbability>::lowest();
                       });
        return likelihoods;
    });
    return result;
}

GenotypeLogMarginalVector
init_genotype_log_marginals(const PopulationModel::GenotypeVector& genotypes,
                            const HardyWeinbergModel& hw_model)
{
    GenotypeLogMarginalVector result {};
    result.reserve(genotypes.size());
    for (const auto& genotype : genotypes) {
        result.push_back({genotype, hw_model.evaluate(genotype)});
    }
    return result;
}

void update_genotype_log_marginals(GenotypeLogMarginalVector& current_log_marginals,
                                   const HardyWeinbergModel& hw_model)
{
    std::for_each(std::begin(current_log_marginals), std::end(current_log_marginals),
                  [&hw_model] (auto& p) { p.log_probability = hw_model.evaluate(p.genotype); });
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
                       std::cbegin(sample_genotype_log_likilhoods), std::begin(posteriors),
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
    auto likelihood_itr = std::cbegin(genotype_log_likilhoods);
    for (auto& sample_genotype_posteriors : current_genotype_posteriors) {
        std::transform(std::cbegin(genotype_log_marginals), std::cend(genotype_log_marginals),
                       std::cbegin(*likelihood_itr++), std::begin(sample_genotype_posteriors),
                       [] (const auto& log_marginal, const auto& log_likeilhood) {
                           return log_marginal.log_probability + log_likeilhood;
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

double update_haplotype_frequencies(HardyWeinbergModel& hw_model,
                                    const GenotypeMarginalPosteriorMatrix& genotype_posteriors,
                                    const InverseGenotypeTable& genotypes_containing_haplotypes,
                                    const std::size_t num_haplotypes,
                                    const double frequency_update_norm)
{
    const auto collaped_posteriors = collapse_genotype_posteriors(genotype_posteriors);
    double max_frequency_change {0};
    auto& current_haplotype_frequencies = hw_model.frequencies();
    for (std::size_t haplotype_idx {0}; haplotype_idx < num_haplotypes; ++haplotype_idx) {
        auto& current_frequency = current_haplotype_frequencies[haplotype_idx];
        double new_frequency {0};
        for (const auto& genotype_index : genotypes_containing_haplotypes[haplotype_idx]) {
            new_frequency += collaped_posteriors[genotype_index];
        }
        new_frequency /= frequency_update_norm;
        const auto frequency_change = std::abs(current_frequency - new_frequency);
        if (frequency_change > max_frequency_change) {
            max_frequency_change = frequency_change;
        }
        current_frequency = std::max(new_frequency, std::numeric_limits<decltype(new_frequency)>::min());
    }
    return max_frequency_change;
}

double do_em_iteration(GenotypeMarginalPosteriorMatrix& genotype_posteriors,
                       HardyWeinbergModel& hw_model,
                       GenotypeLogMarginalVector& genotype_log_marginals,
                       const ModelConstants& constants)
{
    const auto max_change = update_haplotype_frequencies(hw_model,
                                                         genotype_posteriors,
                                                         constants.genotypes_containing_haplotypes,
                                                         constants.num_haplotypes,
                                                         constants.frequency_update_norm);
    update_genotype_log_marginals(genotype_log_marginals, hw_model);
    update_genotype_posteriors(genotype_posteriors, genotype_log_marginals, constants.genotype_log_likilhoods);
    return max_change;
}

void run_em(GenotypeMarginalPosteriorMatrix& genotype_posteriors,
            HardyWeinbergModel& hw_model,
            GenotypeLogMarginalVector& genotype_log_marginals,
            const ModelConstants& constants,
            const EMOptions options,
            boost::optional<logging::TraceLogger> trace_log = boost::none)
{
    for (unsigned n {1}; n <= options.max_iterations; ++n) {
        const auto max_change = do_em_iteration(genotype_posteriors, hw_model, genotype_log_marginals, constants);
        if (max_change <= options.epsilon) break;
    }
}

auto compute_approx_genotype_marginal_posteriors(const PopulationModel::GenotypeVector& genotypes,
                                                 const GenotypeLogLikelihoodMatrix& genotype_likelihoods,
                                                 const ModelConstants& constants,
                                                 const EMOptions options)
{
    auto hw_model = make_hardy_weinberg_model(constants);
    auto genotype_log_marginals = init_genotype_log_marginals(genotypes, hw_model);
    auto result = init_genotype_posteriors(genotype_log_marginals, genotype_likelihoods);
    run_em(result, hw_model, genotype_log_marginals, constants, options);
    return result;
}

auto compute_approx_genotype_marginal_posteriors(const MappableBlock<Haplotype>& haplotypes,
                                                 const PopulationModel::GenotypeVector& genotypes,
                                                 const GenotypeLogLikelihoodMatrix& genotype_likelihoods,
                                                 const EMOptions options)
{
    const ModelConstants constants {haplotypes, genotypes, genotype_likelihoods};
    return compute_approx_genotype_marginal_posteriors(genotypes, genotype_likelihoods, constants, options);
}

auto compute_approx_genotype_marginal_posteriors(const MappableBlock<Haplotype>& haplotypes,
                                                 const PopulationModel::GenotypeVector& genotypes,
                                                 const GenotypeLogLikelihoodMatrix& genotype_likelihoods,
                                                 const std::vector<unsigned>& sample_plodies,
                                                 const EMOptions options)
{
    const ModelConstants constants {haplotypes, genotypes, genotype_likelihoods, sample_plodies};
    return compute_approx_genotype_marginal_posteriors(genotypes, genotype_likelihoods, constants, options);
}

using GenotypeCombinationVector = std::vector<std::size_t>;
using GenotypeCombinationMatrix = std::vector<GenotypeCombinationVector>;

boost::optional<std::size_t> compute_num_combinations(const std::size_t num_genotypes, const std::size_t num_samples)
{
    if (maths::is_safe_ipow(num_genotypes, num_samples)) {
        return maths::ipow(num_genotypes, num_samples);
    } else {
        return boost::none;
    }
}

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

boost::optional<std::size_t> compute_num_combinations(const std::vector<std::size_t>& sample_genotype_set_ids, const std::vector<std::size_t>& genotype_set_sizes)
{
    std::size_t result {1};
    for (std::size_t genotype_set_id {0}; genotype_set_id < genotype_set_sizes.size(); ++genotype_set_id) {
        std::size_t sample_count = std::count(std::cbegin(sample_genotype_set_ids), std::cend(sample_genotype_set_ids), genotype_set_id);
        const auto num_combinations = compute_num_combinations(genotype_set_sizes[genotype_set_id], sample_count);
        if (!num_combinations) {
            return boost::none; // overflow
        }
        const auto new_result = result * *num_combinations;
        if (result != 0 && new_result / result != *num_combinations) {
            return boost::none; // overflow
        }
        result = new_result;
    }
    return result;
}

auto generate_all_genotype_combinations(const std::vector<std::size_t>& sample_genotype_set_ids, const std::vector<std::size_t>& genotype_set_sizes)
{
    const auto num_genotype_sets = genotype_set_sizes.size();
    std::vector<GenotypeCombinationMatrix> combinations {};
    combinations.reserve(num_genotype_sets);
    std::size_t total_num_combinations {1};
    const auto num_samples = sample_genotype_set_ids.size();
    std::vector<std::size_t> sample_mappings(num_samples);
    std::size_t sample_mapping_idx {0};
    for (std::size_t genotype_set_id {0}; genotype_set_id < num_genotype_sets; ++genotype_set_id) {
        std::size_t num_samples_in_genotype_set {0};
        for (std::size_t sample_idx {0}; sample_idx < num_samples; ++sample_idx) {
            if (sample_genotype_set_ids[sample_idx] == genotype_set_id) {
                sample_mappings[sample_idx] = sample_mapping_idx++;
                ++num_samples_in_genotype_set;
            }
        }
        auto set_combinations = generate_all_genotype_combinations(genotype_set_sizes[genotype_set_id], num_samples_in_genotype_set);
        if (genotype_set_id > 0) {
            for (auto& combination : set_combinations) {
                for (auto& genotype_idx : combination) {
                    genotype_idx += genotype_set_sizes[genotype_set_id - 1];
                }
            }
        }
        total_num_combinations *= set_combinations.size();
        combinations.push_back(std::move(set_combinations));
    }
    GenotypeCombinationMatrix result {};
    result.reserve(total_num_combinations);
    for (const auto& set_combinations : combinations) {
        if (!result.empty()) {
            GenotypeCombinationMatrix tmp {};
            tmp.reserve(result.size() * set_combinations.size());
            for (const auto& combination : result) {
                for (const auto& set_combination : set_combinations) {
                    tmp.push_back(concat(combination, set_combination));
                }
            }
            result = std::move(tmp);
        } else {
            result = set_combinations;
        }
    }
    for (auto& combination : result) {
        GenotypeCombinationVector tmp(num_samples);
        for (std::size_t sample_idx {0}; sample_idx < num_samples; ++sample_idx) {
            tmp[sample_idx] = combination[sample_mappings[sample_idx]];
        }
        combination = std::move(tmp);
    }
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

template <typename T>
auto zip_index(const std::vector<T>& v)
{
    std::vector<std::pair<T, unsigned>> result(v.size());
    for (unsigned idx {0}; idx < v.size(); ++idx) {
        result[idx] = std::make_pair(v[idx], idx);
    }
    return result;
}

std::vector<unsigned>
select_top_k_genotypes(const PopulationModel::GenotypeVector& genotypes,
                       const GenotypeMarginalPosteriorMatrix& em_genotype_marginals,
                       const std::size_t k)
{
    if (genotypes.size() <= k) {
        std::vector<unsigned> result(genotypes.size());
        std::iota(std::begin(result), std::end(result), 0);
        return result;
    } else {
        std::vector<std::vector<std::pair<double, unsigned>>> indexed_marginals {};
        indexed_marginals.reserve(em_genotype_marginals.size());
        for (const auto& marginals : em_genotype_marginals) {
            auto tmp = zip_index(marginals);
            std::nth_element(std::begin(tmp), std::next(std::begin(tmp), k), std::end(tmp), std::greater<> {});
            indexed_marginals.push_back(std::move(tmp));
        }
        std::vector<unsigned> result {}, top(genotypes.size(), 0u);
        result.reserve(k);
        for (std::size_t j {0}; j <= k; ++j) {
            for (const auto& marginals : indexed_marginals) {
                ++top[marginals.front().second];
            }
            const auto max_itr = std::max_element(std::begin(top), std::end(top));
            const auto max_idx = static_cast<unsigned>(std::distance(std::begin(top), max_itr));
            if (std::find(std::cbegin(result), std::cend(result), max_idx) == std::cend(result)) {
                result.push_back(max_idx);
            }
            *max_itr = 0;
            for (auto& marginals : indexed_marginals) {
                if (marginals.front().second == max_idx) {
                    marginals.erase(std::cbegin(marginals));
                }
            }
        }
        return result;
    }
}

auto propose_genotype_combinations(const PopulationModel::GenotypeVector& genotypes,
                                   const GenotypeMarginalPosteriorMatrix& em_genotype_marginals,
                                   const std::size_t max_genotype_combinations)
{
    const auto num_samples = em_genotype_marginals.size();
    const auto max_possible_genotype_combinations = compute_num_combinations(genotypes.size(), num_samples);
    if (max_possible_genotype_combinations && *max_possible_genotype_combinations <= max_genotype_combinations) {
        return generate_all_genotype_combinations(genotypes.size(), num_samples);
    }
    auto result = select_top_k_tuples(em_genotype_marginals, max_genotype_combinations);
    const auto top_k_genotype_indices = select_top_k_genotypes(genotypes, em_genotype_marginals, num_samples / 2);
    for (const auto genotype_idx : top_k_genotype_indices) {
        for (std::size_t sample_idx {0}; sample_idx < num_samples; ++sample_idx) {
            if (result.front()[sample_idx] != genotype_idx) {
                auto tmp = result.front();
                tmp[sample_idx] = genotype_idx;
                if (std::find(std::cbegin(result), std::cend(result), tmp) == std::cend(result)) {
                    result.push_back(std::move(tmp));
                }
            }
        }
    }
    const auto hom_ref_idx = find_hom_ref_idx(genotypes);
    if (hom_ref_idx) {
        std::vector<std::size_t> ref_indices(num_samples, *hom_ref_idx);
        if (std::find(std::cbegin(result), std::cend(result), ref_indices) == std::cend(result)) {
            result.back() = std::move(ref_indices);
        }
    }
    return result;
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

template <typename Range, typename V>
void fill(const Range& genotypes, const GenotypeCombinationVector& combination, V& result)
{
    result.clear();
    std::transform(std::cbegin(combination), std::cend(combination), std::back_inserter(result),
                   [&] (const auto index) { return std::cref(genotypes[index]); });
}

auto calculate_posteriors(const PopulationModel::GenotypeVector& genotypes,
                          const GenotypeCombinationMatrix& genotype_combinations,
                          const GenotypeLogLikelihoodMatrix& genotype_likelihoods,
                          const PopulationPriorModel& prior_model)
{
    std::vector<double> result {};
    GenotypeLogLikelihoodVector likelihoods_buffer(genotype_likelihoods.size());
    PopulationPriorModel::GenotypeReferenceVector genotype_combination {};
    genotype_combination.reserve(genotype_combinations.front().size());
    for (const auto& combination : genotype_combinations) {
        fill(genotype_likelihoods, combination, likelihoods_buffer);
        fill(genotypes, combination, genotype_combination);
        result.push_back(prior_model.evaluate(genotype_combination) + sum(likelihoods_buffer));
    }
    const auto norm = maths::normalise_exp(result);
    return std::make_pair(std::move(result), norm);
}

void set_posterior_marginals(const GenotypeCombinationMatrix& genotype_combinations,
                             const std::vector<double>& joint_posteriors,
                             const std::size_t num_genotypes, const std::size_t num_samples,
                             PopulationModel::InferredLatents& result)
{
    assert(joint_posteriors.size() == genotype_combinations.size());
    std::vector<std::vector<double>> marginals(num_samples, std::vector<double>(num_genotypes, 0.0));
    for (std::size_t i {0}; i < genotype_combinations.size(); ++i) {
        assert(genotype_combinations[i].size() == num_samples);
        for (std::size_t s {0}; s < num_samples; ++s) {
            marginals[s][genotype_combinations[i][s]] += joint_posteriors[i];
        }
    }
    result.posteriors.marginal_genotype_probabilities = std::move(marginals);
}

template <typename Range>
void calculate_posterior_marginals(const Range& genotypes,
                                   const GenotypeCombinationMatrix& genotype_combinations,
                                   const GenotypeLogLikelihoodMatrix& genotype_likelihoods,
                                   const PopulationPriorModel& prior_model,
                                   PopulationModel::InferredLatents& result)
{
    std::vector<double> joint_posteriors; double norm;
    std::tie(joint_posteriors, norm) = calculate_posteriors(genotypes, genotype_combinations, genotype_likelihoods, prior_model);
    const auto num_samples = genotype_likelihoods.size();
    set_posterior_marginals(genotype_combinations, joint_posteriors, genotypes.size(), num_samples, result);
    result.log_evidence = norm;
}

} // namespace

PopulationModel::InferredLatents
PopulationModel::evaluate(const SampleVector& samples,
                          const MappableBlock<Haplotype>& haplotypes,
                          const GenotypeVector& genotypes,
                          const HaplotypeLikelihoodArray& haplotype_likelihoods) const
{
    assert(!genotypes.empty());
    const auto genotype_log_likelihoods = compute_genotype_log_likelihoods(samples, genotypes, haplotype_likelihoods);
    const auto num_possible_genotype_combinations = compute_num_combinations(genotypes.size(), samples.size());
    InferredLatents result;
    GenotypeCombinationMatrix genotype_combinations {};
    if (!options_.max_genotype_combinations || (num_possible_genotype_combinations && *num_possible_genotype_combinations <= *options_.max_genotype_combinations)) {
        genotype_combinations = generate_all_genotype_combinations(genotypes.size(), samples.size());
    } else {
        const auto max_genotype_combinations = options_.max_genotype_combinations ? *options_.max_genotype_combinations : *num_possible_genotype_combinations;
        const EMOptions em_options {options_.max_em_iterations, options_.em_epsilon};
        const auto em_genotype_marginals = compute_approx_genotype_marginal_posteriors(haplotypes, genotypes, genotype_log_likelihoods, em_options);
        genotype_combinations = propose_genotype_combinations(genotypes, em_genotype_marginals, max_genotype_combinations);
    }
    calculate_posterior_marginals(genotypes, genotype_combinations, genotype_log_likelihoods, prior_model_, result);
    return result;
}

namespace {

using GenotypeMask = std::vector<bool>;

GenotypeMask
make_genotype_mask(const unsigned ploidy, const PopulationModel::GenotypeVector& genotypes)
{
    GenotypeMask result(genotypes.size());
    const auto is_ploidy = [ploidy] (const auto& genotype) { return genotype.ploidy() == ploidy; };
    std::transform(std::cbegin(genotypes), std::cend(genotypes), std::begin(result), is_ploidy);
    return result;
}

std::vector<GenotypeMask>
make_genotype_masks(const std::vector<unsigned>& sample_ploidies,
                    const PopulationModel::GenotypeVector& genotypes)
{
    std::vector<GenotypeMask> result {};
    result.reserve(sample_ploidies.size());
    for (auto ploidy : sample_ploidies) result.push_back(make_genotype_mask(ploidy, genotypes));
    return result;
}

unsigned find_max_ploidy(const PopulationModel::GenotypeVector& genotypes) noexcept
{
    const static auto ploidy_less = [] (const auto& lhs, const auto& rhs) { return lhs.ploidy() < rhs.ploidy(); };
    return std::max_element(std::cbegin(genotypes), std::cend(genotypes), ploidy_less)->ploidy();
}

std::pair<std::vector<std::size_t>, std::vector<std::size_t>>
get_genotype_sets(const std::vector<unsigned>& sample_ploidies,
                  const PopulationModel::GenotypeVector& genotypes)
{
    const auto max_ploidy = find_max_ploidy(genotypes);
    std::vector<std::size_t> genotype_ploidy_sizes(max_ploidy + 1);
    for (const auto& genotype : genotypes) {
        ++genotype_ploidy_sizes[genotype.ploidy()];
    }
    const auto num_genotype_sets = std::count_if(std::cbegin(genotype_ploidy_sizes), std::cend(genotype_ploidy_sizes), [] (auto n) { return n > 0; });
    std::vector<std::size_t> genotype_set_sizes(num_genotype_sets), ploidy_to_genotype_set_idx(max_ploidy + 1);
    for (unsigned ploidy {0}, genotype_set_idx {0}; ploidy <= max_ploidy; ++ploidy) {
        if (genotype_ploidy_sizes[ploidy] > 0) {
            ploidy_to_genotype_set_idx[ploidy] = genotype_set_idx;
            genotype_set_sizes[genotype_set_idx++] = genotype_ploidy_sizes[ploidy];
        }
    }
    std::vector<std::size_t> sample_genotype_set_ids(sample_ploidies.size());
    for (std::size_t sample_idx {0}; sample_idx < sample_ploidies.size(); ++sample_idx) {
        sample_genotype_set_ids[sample_idx] = ploidy_to_genotype_set_idx[sample_ploidies[sample_idx]];
    }
    return {std::move(sample_genotype_set_ids), std::move(genotype_set_sizes)};
}

} // namespace

PopulationModel::InferredLatents
PopulationModel::evaluate(const SampleVector& samples,
                          const std::vector<unsigned>& sample_ploidies,
                          const MappableBlock<Haplotype>& haplotypes,
                          const GenotypeVector& genotypes,
                          const HaplotypeLikelihoodArray& haplotype_likelihoods) const
{
    const auto genotype_masks = make_genotype_masks(sample_ploidies, genotypes);
    const auto genotype_log_likelihoods = compute_genotype_log_likelihoods(samples, genotypes, haplotype_likelihoods, genotype_masks);
    std::vector<std::size_t> sample_genotype_set_ids, genotype_set_sizes;
    std::tie(sample_genotype_set_ids, genotype_set_sizes) = get_genotype_sets(sample_ploidies, genotypes);
    const auto num_possible_genotype_combinations = compute_num_combinations(sample_genotype_set_ids, genotype_set_sizes);
    InferredLatents result {};
    GenotypeCombinationMatrix genotype_combinations {};
    if (!options_.max_genotype_combinations || (num_possible_genotype_combinations && *num_possible_genotype_combinations <= *options_.max_genotype_combinations)) {
        genotype_combinations = generate_all_genotype_combinations(sample_genotype_set_ids, genotype_set_sizes);
    } else {
        const auto max_genotype_combinations = options_.max_genotype_combinations ? *options_.max_genotype_combinations : *num_possible_genotype_combinations;
        const EMOptions em_options {options_.max_em_iterations, options_.em_epsilon};
        const auto em_genotype_marginals = compute_approx_genotype_marginal_posteriors(haplotypes, genotypes, genotype_log_likelihoods, sample_ploidies, em_options);
        genotype_combinations = propose_genotype_combinations(genotypes, em_genotype_marginals, max_genotype_combinations);
    }
    calculate_posterior_marginals(genotypes, genotype_combinations, genotype_log_likelihoods, prior_model_, result);
    return result;
}

namespace debug {
    
} // namespace debug

} // namespace model
} // namespace octopus
