// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "population_model.hpp"

#include <unordered_map>
#include <utility>
#include <algorithm>
#include <cmath>
#include <cassert>
#include <iostream>

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
        auto itr = result_map.emplace(std::piecewise_construct,
                                      std::forward_as_tuple(std::cref(haplotype)),
                                      std::forward_as_tuple());
        itr.first->second.reserve(cardinality);
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

using HaplotypeReference    = std::reference_wrapper<const Haplotype>;
using HaplotypeFrequencyMap = std::unordered_map<HaplotypeReference, double>;

double calculate_frequency_update_norm(const std::size_t num_samples, const unsigned ploidy)
{
    return static_cast<double>(num_samples) * ploidy;
}

struct EmOptions
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
            const ModelConstants& constants, const EmOptions options,
            boost::optional<logging::TraceLogger> trace_log = boost::none)
{
    for (unsigned n {1}; n <= options.max_iterations; ++n) {
        const auto max_change = do_em_iteration(genotype_posteriors, haplotype_frequencies,
                                                genotype_log_marginals,constants);
        if (max_change <= options.epsilon) break;
    }
}

auto compute_approx_genotype_marginal_posteriors(const std::vector<Genotype<Haplotype>>& genotypes,
                                                 const GenotypeLogLikelihoodMatrix& genotype_likelihoods,
                                                 const EmOptions options)
{
    const auto haplotypes = extract_unique_elements(genotypes);
    const ModelConstants constants {haplotypes, genotypes, genotype_likelihoods};
    auto haplotype_frequencies = init_haplotype_frequencies(constants);
    auto genotype_log_marginals = init_genotype_log_marginals(genotypes, haplotype_frequencies);
    auto result = init_genotype_posteriors(genotype_log_marginals, genotype_likelihoods);
    run_em(result, haplotype_frequencies, genotype_log_marginals, constants, options);
    return result;
}

using GenotypeCombinationVector = std::vector<std::size_t>;
using GenotypeCombinationMatrix = std::vector<GenotypeCombinationVector>;

auto num_combinations(const std::size_t num_genotypes, const std::size_t num_samples)
{
    return std::pow(num_genotypes, num_samples);
}

auto get_all_genotype_combinations(const std::size_t num_genotypes, const std::size_t num_samples)
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

auto index_and_sort(const GenotypeMarginalPosteriorVector& genotype_posteriors, const std::size_t k)
{
    auto result = index(genotype_posteriors);
    const auto middle = std::next(std::begin(result), std::min(k, result.size()));
    std::partial_sort(std::begin(result), middle, std::end(result), std::greater<> {});
    std::for_each(std::begin(result), middle, [] (auto& p) { p.first = std::log(p.first); });
    result.erase(middle, std::end(result));
    return result;
}

using IndexedProbability = std::pair<double, std::size_t>;
using IndexedProbabilityVector = std::vector<IndexedProbability>;

struct CombinationProbabilityRow
{
    std::vector<std::size_t> combination;
    double log_probability;
};

using CombinationProbabilityMatrix = std::vector<CombinationProbabilityRow>;

auto get_differences(const IndexedProbabilityVector& genotype_posteriors)
{
    std::vector<double> result(genotype_posteriors.size());
    std::transform(std::cbegin(genotype_posteriors), std::cend(genotype_posteriors), std::begin(result),
                   [] (const auto& p) { return p.first; });
    std::adjacent_difference(std::cbegin(result), std::cend(result), std::begin(result));
    return result;
}

auto get_differences(const CombinationProbabilityMatrix& matrix)
{
    std::vector<double> result(matrix.size());
    std::transform(std::cbegin(matrix), std::cend(matrix), std::begin(result),
                   [] (const auto& p) { return p.log_probability; });
    std::adjacent_difference(std::cbegin(result), std::cend(result), std::begin(result));
    return result;
}

void join(const IndexedProbabilityVector& genotype_posteriors,
          CombinationProbabilityMatrix& result,
          const std::size_t k)
{
    const auto n = std::min(k, genotype_posteriors.size());
    if (result.empty()) {
        std::transform(std::cbegin(genotype_posteriors), std::next(std::cbegin(genotype_posteriors), n),
                       std::back_inserter(result), [=] (const auto& p) -> CombinationProbabilityRow {
                           return {{p.second}, p.first};
                       });
    } else {
        const auto m = result.size();
//        const auto differences1 = get_differences(result);
//        const auto differences2 = get_differences(genotype_posteriors);
        const auto K = std::min(k, n * m);
        CombinationProbabilityMatrix tmp {};
        tmp.reserve(n * m);
        for (std::size_t i {0}; i < m; ++i) {
            for (std::size_t j {0}; j < n; ++j) {
                tmp.push_back(result[i]);
                tmp.back().combination.push_back(genotype_posteriors[j].second);
                tmp.back().log_probability += genotype_posteriors[j].first;
            }
        }
        const auto middle = std::next(std::begin(tmp), K);
        std::partial_sort(std::begin(tmp), middle, std::end(tmp),
                          [] (const auto& lhs, const auto& rhs) {
                              return lhs.log_probability > rhs.log_probability;
                          });
        tmp.erase(middle, std::end(tmp));
//        tmp.reserve(K);
//        for (std::size_t i {0}, j {0}, t {0}; t < K; ++t) {
//            assert(i < m && j < n);
//            std::cout << i << " " << j << std::endl;
//            tmp.push_back(result[i]);
//            tmp.back().combination.push_back(genotype_posteriors[j].second);
//            tmp.back().log_probability += genotype_posteriors[j].first;
//            const auto prev_score = tmp.back().log_probability;
//            if (j < n - 1 && (i == m - 1 || differences1[i + 1] < differences2[j + 1])) {
//                ++j;
//                while (i > 0 && result[i - 1].log_probability + genotype_posteriors[j].first < prev_score) {
//                    --i;
//                }
//            } else {
//                ++i;
//                while (j > 0 && result[i].log_probability + genotype_posteriors[j - 1].first < prev_score) {
//                    --j;
//                }
//            }
//        }
        result = std::move(tmp);
    }
}

auto get_genotype_combinations(const std::vector<Genotype<Haplotype>>& genotypes,
                               const GenotypeMarginalPosteriorMatrix& genotype_posteriors,
                               const std::size_t max_combinations)
{
    const auto num_samples = genotype_posteriors.size();
    assert(max_combinations >= num_samples);
    const auto num_possible_combinations = num_combinations(genotypes.size(), num_samples);
    if (num_possible_combinations <= max_combinations) {
        return get_all_genotype_combinations(genotypes.size(), num_samples);
    }
    CombinationProbabilityMatrix combinations {};
    combinations.reserve(max_combinations);
    for (const auto& sample : genotype_posteriors) {
        join(index_and_sort(sample, max_combinations), combinations, max_combinations);
    }
    GenotypeCombinationMatrix result {};
    result.reserve(max_combinations);
    for (auto&& row : combinations) {
        result.push_back(std::move(row.combination));
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
                          const GenotypeCombinationMatrix& genotype_combinations,
                          const GenotypeLogLikelihoodMatrix& genotype_likelihoods,
                          const PopulationPriorModel& prior_model)
{
    assert(!genotypes.empty());
    std::vector<double> result {};
    GenotypeLogLikelihoodVector tmp_likelihoods(genotype_likelihoods.size());
    GenotypeReferenceVector tmp_genotypes {};
    for (const auto& indices : genotype_combinations) {
        fill(genotype_likelihoods, indices, tmp_likelihoods);
        fill(genotypes, indices, tmp_genotypes);
        result.push_back(prior_model.evaluate(tmp_genotypes) + sum(tmp_likelihoods));
    }
    const auto norm = maths::normalise_exp(result);
    return std::make_pair(result, norm);
}

} // namespace

PopulationModel::InferredLatents
PopulationModel::evaluate(const SampleVector& samples, const GenotypeVector& genotypes,
                          const HaplotypeLikelihoodCache& haplotype_likelihoods) const
{
    assert(!genotypes.empty());
    const auto genotype_log_likelihoods = compute_genotype_log_likelihoods(samples, genotypes, haplotype_likelihoods);
    const auto approx_genotype_posteriors = compute_approx_genotype_marginal_posteriors(genotypes, genotype_log_likelihoods,
                                                                                        {options_.max_em_iterations, 0.0001});
    const auto max_combinations = options_.max_combinations_per_sample * samples.size();
    auto genotype_combinations = get_genotype_combinations(genotypes, approx_genotype_posteriors, max_combinations);
    auto p = calculate_posteriors(genotypes, genotype_combinations, genotype_log_likelihoods, prior_model_);
    return {{std::move(genotype_combinations), std::move(p.first)}, p.second};
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
