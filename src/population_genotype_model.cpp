//
//  population_genotype_model.cpp
//  Octopus
//
//  Created by Daniel Cooke on 26/08/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "population_genotype_model.hpp"

#include <algorithm>
#include <numeric>
#include <cmath>
#include <cassert>
#include <iostream>

#include "hardy_weinberg_model.hpp"
#include "dirichlet_model.hpp"
#include "fixed_ploidy_genotype_likelihood_model.hpp"
#include "maths.hpp"
#include "read_utils.hpp"
#include "logging.hpp"

namespace Octopus
{
namespace GenotypeModel
{
// public methods

Population::Population(const unsigned ploidy, const unsigned max_em_iterations, const double em_epsilon)
:
ploidy_ {ploidy},
max_em_iterations_ {max_em_iterations},
em_epsilon_ {em_epsilon}
{}

// non member methods

namespace
{
using GenotypeLogLikelihood        = double;
using SampleGenotypeLogLikelihoods = std::vector<GenotypeLogLikelihood>;
using GenotypeLogLikelihoodMap     = std::unordered_map<SampleIdType, SampleGenotypeLogLikelihoods>;

struct GenotypeLogProbability
{
    GenotypeLogProbability() = delete;
    GenotypeLogProbability(const Genotype<Haplotype>& genotype, double log_probability)
    : genotype {genotype}, log_probability {log_probability} {}
    
    const Genotype<Haplotype>& genotype;
    double log_probability;
};

using GenotypeLogMarginals = std::vector<GenotypeLogProbability>;

using GenotypePosterior        = double;
using SampleGenotypePosteriors = std::vector<GenotypePosterior>;
using GenotypePosteriorMap     = std::unordered_map<SampleIdType, SampleGenotypePosteriors>;

using HaplotypePriorCounts = std::vector<double>;

HaplotypePriorCounts flatten_haplotype_prior_counts(const std::vector<Haplotype>& haplotypes,
                                                    const HaplotypePriorCountMap& prior_counts)
{
    HaplotypePriorCounts result(haplotypes.size());
    
    std::transform(std::cbegin(haplotypes), std::cend(haplotypes), std::begin(result),
                   [&prior_counts] (const auto& haplotype) {
                       return prior_counts.at(haplotype);
                   });
    
    return result;
}

auto make_inverse_genotype_table(const std::vector<Haplotype>& haplotypes,
                                 const std::vector<Genotype<Haplotype>>& genotypes)
{
    assert(!haplotypes.empty() && !genotypes.empty());
    
    using HaplotypeReference = std::reference_wrapper<const Haplotype>;
    
    std::unordered_map<HaplotypeReference, std::vector<std::size_t>> result_map {haplotypes.size()};
    
    const auto r = element_cardinality_in_genotypes(static_cast<unsigned>(haplotypes.size()),
                                                    genotypes.front().ploidy());
    
    for (const auto& haplotype : haplotypes) {
        auto it = result_map.emplace(std::piecewise_construct,
                                     std::forward_as_tuple(std::cref(haplotype)),
                                     std::forward_as_tuple());
        it.first->second.reserve(r);
    }
    
    for (std::size_t i {0}; i < genotypes.size(); ++i) {
        for (const auto& haplotype : genotypes[i]) {
            result_map.at(haplotype).emplace_back(i);
        }
    }
    
    std::vector<std::vector<std::size_t>> result {};
    
    for (const auto& haplotype : haplotypes) {
        result.emplace_back(std::move(result_map.at(haplotype)));
    }
    
    return result;
}

double calculate_frequency_update_norm(const std::size_t num_samples, const unsigned ploidy,
                                       const std::size_t num_haplotypes, const double& prior_count_sum)
{
    return static_cast<double>(num_samples) * ploidy + prior_count_sum;
               // - static_cast<double>(num_haplotypes);
}
    
struct ModelConstants
{
    using InverseGenotypeTable = std::vector<std::vector<std::size_t>>;
    
    const std::vector<Haplotype>& haplotypes;
    const std::vector<Genotype<Haplotype>>& genotypes;
    
    const GenotypeLogLikelihoodMap& genotype_log_likilhoods;
    const HaplotypePriorCounts haplotype_prior_counts;
    
    const unsigned ploidy;
    
    const double prior_count_sum;
    const double frequency_update_norm;
    
    const InverseGenotypeTable genotypes_containing_haplotypes;
    
    ModelConstants() = delete;
    
    explicit ModelConstants(const std::vector<Haplotype>& haplotypes,
                            const std::vector<Genotype<Haplotype>>& genotypes,
                            const GenotypeLogLikelihoodMap& genotype_log_likilhoods,
                            const HaplotypePriorCountMap& haplotype_prior_counts)
    :
    haplotypes {haplotypes},
    genotypes {genotypes},
    genotype_log_likilhoods {genotype_log_likilhoods},
    haplotype_prior_counts {flatten_haplotype_prior_counts(haplotypes, haplotype_prior_counts)},
    ploidy {genotypes.front().ploidy()},
    prior_count_sum {std::accumulate(std::cbegin(this->haplotype_prior_counts),
                                     std::cend(this->haplotype_prior_counts), 0.0)},
    frequency_update_norm {calculate_frequency_update_norm(genotype_log_likilhoods.size(), ploidy,
                                                           haplotypes.size(), prior_count_sum)},
    genotypes_containing_haplotypes {make_inverse_genotype_table(haplotypes, genotypes)}
    {
        if (TRACE_MODE) {
            Logging::TraceLogger log {};
            stream(log) << "Prior count sum = " << prior_count_sum;
            stream(log) << "Frequency update norm = " << frequency_update_norm;
        }
    }
    
    ~ModelConstants() = default;
};
} // namespace

namespace debug
{
    void print_genotypes(const std::vector<Genotype<Haplotype>>& genotypes);
    void print_haplotype_priors(const HaplotypePriorCountMap& prior_counts,
                                std::size_t n = 5);
    template <typename S>
    void print_genotype_log_likelihoods(S&&stream,
                                        const std::vector<Genotype<Haplotype>>& genotypes,
                                        const GenotypeLogLikelihoodMap& log_likelihoods,
                                        std::size_t n = 5);
    void print_genotype_log_likelihoods(const std::vector<Genotype<Haplotype>>& genotypes,
                                        const GenotypeLogLikelihoodMap& log_likelihoods,
                                        std::size_t n = 5);
    template <typename S>
    void print_haplotype_frequencies(S&& stream,
                                     const HaplotypeFrequencyMap& haplotype_frequencies,
                                     const std::size_t n = 5);
    void print_haplotype_frequencies(const HaplotypeFrequencyMap& haplotype_frequencies,
                                     std::size_t n = 5);
    template <typename S>
    void print_genotype_log_marginals(S&& stream,
                                      const std::vector<Genotype<Haplotype>>& genotypes,
                                      const GenotypeLogMarginals& genotype_log_marginals,
                                      std::size_t n = 5);
    void print_genotype_log_marginals(const std::vector<Genotype<Haplotype>>& genotypes,
                                      const GenotypeLogMarginals& genotype_log_marginals,
                                      std::size_t n = 5);
    template <typename S>
    void print_genotype_posteriors(S&& stream,
                                   const std::vector<Genotype<Haplotype>>& genotypes,
                                   const GenotypePosteriorMap& genotype_posteriors,
                                   std::size_t n = 5);
    void print_genotype_posteriors(const std::vector<Genotype<Haplotype>>& genotypes,
                                   const GenotypePosteriorMap& genotype_posteriors,
                                   std::size_t n = 5);
} // namespace debug

// non member methods

namespace
{
HaplotypeFrequencyMap
init_haplotype_frequencies(const HaplotypePriorCountMap& haplotype_prior_counts,
                           const ModelConstants& constants)
{
    HaplotypeFrequencyMap result {haplotype_prior_counts.size()};
    
    for (const auto& haplotype_count : haplotype_prior_counts) {
        result.emplace(haplotype_count.first, haplotype_count.second / constants.prior_count_sum);
    }
    
    return result;
}

GenotypeLogLikelihoodMap
compute_genotype_log_likelihoods(const std::vector<SampleIdType>& samples,
                                 const std::vector<Genotype<Haplotype>>& genotypes,
                                 const HaplotypeLikelihoodCache& haplotype_likelihoods)
{
    assert(!genotypes.empty());
    
    FixedPloidyGenotypeLikelihoodModel likelihood_model {genotypes.front().ploidy(), haplotype_likelihoods};
    
    GenotypeLogLikelihoodMap result {samples.size()};
    
    for (const auto& sample : samples) {
        const auto it = result.emplace(std::piecewise_construct, std::forward_as_tuple(sample),
                                       std::forward_as_tuple(genotypes.size())).first;
        
        std::transform(std::cbegin(genotypes), std::cend(genotypes), std::begin(it->second),
                       [&sample, &likelihood_model] (const auto& genotype) {
                           return likelihood_model.log_likelihood(sample, genotype);
                       });
    }
    
    return result;
}

GenotypeLogMarginals
init_genotype_log_marginals(const std::vector<Genotype<Haplotype>>& genotypes,
                            const HaplotypeFrequencyMap& haplotype_frequencies)
{
    GenotypeLogMarginals result {};
    result.reserve(genotypes.size());
    
    for (const auto& genotype : genotypes) {
        result.emplace_back(genotype, log_hardy_weinberg(genotype, haplotype_frequencies));
    }
    
    return result;
}

void update_genotype_log_marginals(GenotypeLogMarginals& current_log_marginals,
                                   const HaplotypeFrequencyMap& haplotype_frequencies)
{
    std::for_each(std::begin(current_log_marginals), std::end(current_log_marginals),
                  [&haplotype_frequencies] (auto& p) {
                      p.log_probability = log_hardy_weinberg(p.genotype, haplotype_frequencies);
                  });
}

void normalise_exp(SampleGenotypePosteriors& unnormalised_log_probabilities)
{
    const auto norm = Maths::log_sum_exp(unnormalised_log_probabilities);
    
    std::transform(std::cbegin(unnormalised_log_probabilities),
                   std::cend(unnormalised_log_probabilities),
                   std::begin(unnormalised_log_probabilities),
                   [norm] (const double p) {
                       return std::exp(p - norm);
                   });
}

GenotypePosteriorMap
init_genotype_posteriors(const GenotypeLogMarginals& genotype_log_marginals,
                         const GenotypeLogLikelihoodMap& genotype_log_likilhoods)
{
    GenotypePosteriorMap result {genotype_log_likilhoods.size()};
    
    for (const auto& sample_genotype_log_likilhoods : genotype_log_likilhoods) {
        SampleGenotypePosteriors sample_result(genotype_log_marginals.size());
        
        std::transform(std::cbegin(genotype_log_marginals), std::cend(genotype_log_marginals),
                       std::cbegin(sample_genotype_log_likilhoods.second),
                       std::begin(sample_result),
                       [] (const auto& genotype_log_marginal, const auto genotype_log_likilhood) {
                           return genotype_log_marginal.log_probability + genotype_log_likilhood;
                       });
        
        normalise_exp(sample_result);
        
        result.emplace(sample_genotype_log_likilhoods.first, std::move(sample_result));
    }
    
    return result;
}

void update_genotype_posteriors(GenotypePosteriorMap& current_genotype_posteriors,
                                const GenotypeLogMarginals& genotype_log_marginals,
                                const GenotypeLogLikelihoodMap& genotype_log_likilhoods)
{
    for (auto& sample_genotype_posteriors : current_genotype_posteriors) {
        std::transform(std::cbegin(genotype_log_marginals), std::cend(genotype_log_marginals),
                       std::cbegin(genotype_log_likilhoods.at(sample_genotype_posteriors.first)),
                       std::begin(sample_genotype_posteriors.second),
                       [] (const auto& log_marginal, const auto& log_likilhood) {
                           return log_marginal.log_probability + log_likilhood;
                       });
        normalise_exp(sample_genotype_posteriors.second);
    }
}

auto collapse_genotype_posteriors(const GenotypePosteriorMap& genotype_posteriors)
{
    assert(!genotype_posteriors.empty());
    
    if (genotype_posteriors.size() == 1) {
        return std::cbegin(genotype_posteriors)->second;
    }
    
    std::vector<double> result(std::cbegin(genotype_posteriors)->second.size());
    
    for (const auto& sample_posteriors : genotype_posteriors) {
        std::transform(std::cbegin(result), std::cend(result),
                       std::cbegin(sample_posteriors.second),
                       std::begin(result),
                       [] (const auto curr, const auto p) {
                           return curr + p;
                       });
    }
    
    return result;
}

double update_haplotype_frequencies(const std::vector<Haplotype>& haplotypes,
                                    HaplotypeFrequencyMap& current_haplotype_frequencies,
                                    const HaplotypePriorCounts& haplotype_prior_counts,
                                    const GenotypePosteriorMap& genotype_posteriors,
                                    const ModelConstants::InverseGenotypeTable& genotypes_containing_haplotypes,
                                    const double frequency_update_norm)
{
    const auto collaped_posteriors = collapse_genotype_posteriors(genotype_posteriors);
    
    double max_frequency_change {0};
    
    for (std::size_t i {0}; i < haplotypes.size(); ++i) {
        auto& current_frequency = current_haplotype_frequencies.at(haplotypes[i]);
        
        double new_frequency {haplotype_prior_counts[i]};
        
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

double do_em_iteration(GenotypePosteriorMap& genotype_posteriors,
                       HaplotypeFrequencyMap& haplotype_frequencies,
                       GenotypeLogMarginals& genotype_log_marginals,
                       const ModelConstants& constants)
{
    const auto max_change = update_haplotype_frequencies(constants.haplotypes,
                                                         haplotype_frequencies,
                                                         constants.haplotype_prior_counts,
                                                         genotype_posteriors,
                                                         constants.genotypes_containing_haplotypes,
                                                         constants.frequency_update_norm);
    
    update_genotype_log_marginals(genotype_log_marginals, haplotype_frequencies);
    
    update_genotype_posteriors(genotype_posteriors, genotype_log_marginals,
                               constants.genotype_log_likilhoods);
    
    return max_change;
}

void run_em(GenotypePosteriorMap& genotype_posteriors,
            HaplotypeFrequencyMap& haplotype_frequencies,
            GenotypeLogMarginals& genotype_log_marginals,
            const ModelConstants& constants,
            const unsigned max_iterations, const double epsilon)
{
    if (TRACE_MODE) {
        Logging::TraceLogger log {};
        stream(log) << "Running EM";
        debug::print_haplotype_frequencies(stream(log), haplotype_frequencies, -1);
        debug::print_genotype_log_marginals(stream(log), constants.genotypes, genotype_log_marginals, -1);
//        debug::print_genotype_posteriors(stream(log), constants.genotypes, genotype_posteriors, -1);
    }
    
    for (unsigned n {1}; n <= max_iterations; ++n) {
        const auto max_change = do_em_iteration(genotype_posteriors, haplotype_frequencies,
                                                genotype_log_marginals,constants);
        
        if (TRACE_MODE) {
            Logging::TraceLogger log {};
            stream(log) << "EM iteration " << n << " : max change = " << max_change;
            debug::print_haplotype_frequencies(stream(log), haplotype_frequencies, -1);
//            debug::print_genotype_log_marginals(stream(log), constants.genotypes, genotype_log_marginals, -1);
//            debug::print_genotype_posteriors(stream(log), constants.genotypes, genotype_posteriors, -1);
        }
        
        if (max_change <= epsilon) break;
    }
}

Population::Latents
make_single_genotype_latents(const std::vector<SampleIdType>& samples,
                             const Genotype<Haplotype>& genotype)
{
    Population::Latents::GenotypeProbabilityMap result {};
    result.push_back(genotype);
    
    for (const auto& sample : samples) {
        insert_sample(sample, std::vector<double> {1}, result);
    }
    
    Population::Latents::HaplotypeFrequencyMap haplotype_frequencies {{genotype[0], 1.0}};
    Population::Latents::HaplotypeProbabilityMap haplotype_posteriors {{genotype[0], 1.0}};
    
    return Population::Latents {
        std::move(haplotype_frequencies),
        std::move(result),
        std::move(haplotype_posteriors)
    };
}

auto calculate_haplotype_posteriors(const std::vector<Haplotype>& haplotypes,
                                    const std::vector<Genotype<Haplotype>>& genotypes,
                                    const SampleGenotypePosteriors& genotype_posteriors)
{
    std::unordered_map<std::reference_wrapper<const Haplotype>, double> result {haplotypes.size()};
    
    for (const auto& haplotype : haplotypes) {
        result.emplace(std::cref(haplotype), 0.0);
    }
    
    auto it = std::cbegin(genotype_posteriors);
    
    for (const auto& genotype : genotypes) {
        for (const auto& haplotype : genotype.copy_unique_ref()) {
            result.at(haplotype) += *it;
        }
        ++it;
    }
    
    return result;
}

auto calculate_haplotype_posteriors(const std::vector<Haplotype>& haplotypes,
                                    const std::vector<Genotype<Haplotype>>& genotypes,
                                    const GenotypePosteriorMap& genotype_posteriors,
                                    const ModelConstants::InverseGenotypeTable& inverse_genotypes)
{
    if (genotype_posteriors.size() == 1) {
        // can optimise this case
        return calculate_haplotype_posteriors(haplotypes, genotypes,
                                              std::cbegin(genotype_posteriors)->second);
    }
    
    std::unordered_map<std::reference_wrapper<const Haplotype>, double> result {haplotypes.size()};
    
    auto it = std::cbegin(inverse_genotypes);
    
    std::vector<std::size_t> genotype_indices(genotypes.size());
    
    std::iota(std::begin(genotype_indices), std::end(genotype_indices), 0);
    
    // noncontaining genotypes are genotypes that do not contain a particular haplotype.
    const auto num_noncontaining_genotypes = genotypes.size() - inverse_genotypes.front().size();
    
    std::vector<std::size_t> noncontaining_genotype_indices(num_noncontaining_genotypes);
    
    for (const auto& haplotype : haplotypes) {
        std::set_difference(std::cbegin(genotype_indices), std::cend(genotype_indices),
                            std::cbegin(*it), std::cend(*it),
                            std::begin(noncontaining_genotype_indices));
        
        double prob_not_observed {1};
        
        for (const auto& sample_genotype_posteriors : genotype_posteriors) {
            prob_not_observed *= std::accumulate(std::cbegin(noncontaining_genotype_indices),
                                                 std::cend(noncontaining_genotype_indices),
                                                 0.0, [&sample_genotype_posteriors]
                                                 (const auto curr, const auto i) {
                                                     return curr + sample_genotype_posteriors.second[i];
                                                 });
        }
        
        result.emplace(std::cref(haplotype), 1.0 - prob_not_observed);
        
        ++it;
    }
    
    return result;
}

Population::Latents
make_latents(const std::vector<Haplotype>& haplotypes,
         std::vector<Genotype<Haplotype>>&& genotypes,
         GenotypePosteriorMap&& genotype_posteriors,
         HaplotypeFrequencyMap&& haplotype_frequencies,
         const ModelConstants& constants)
{
    auto haplotype_posteriors = calculate_haplotype_posteriors(haplotypes, genotypes,
                                                               genotype_posteriors,
                                                               constants.genotypes_containing_haplotypes);
    
    Population::Latents::GenotypeProbabilityMap result {
        std::make_move_iterator(std::begin(genotypes)),
        std::make_move_iterator(std::end(genotypes))
    };
    
    insert_samples(std::move(genotype_posteriors), result);
    
    return Population::Latents {
        std::move(haplotype_frequencies),
        std::move(result),
        std::move(haplotype_posteriors)
    };
}
} // namespace

// private methods

Population::Latents
Population::infer_latents(const std::vector<SampleIdType>& samples,
                          const std::vector<Haplotype>& haplotypes,
                          const HaplotypePrioMap& haplotype_priors,
                          const HaplotypeLikelihoodCache& haplotype_likelihoods) const
{
    assert(!haplotypes.empty());
    
    Logging::DebugLogger log {};
    
    auto genotypes = generate_all_genotypes(haplotypes, ploidy_);
    
    if (DEBUG_MODE) {
        stream(log) << "There are " << genotypes.size() << " genotypes";
    }
    
    assert(!genotypes.empty());
    
    if (genotypes.size() == 1) {
        return make_single_genotype_latents(samples, genotypes.front());
    }
    
    const auto genotype_log_likilhoods = compute_genotype_log_likelihoods(samples, genotypes,
                                                                          haplotype_likelihoods);
    
    if (TRACE_MODE) {
        Logging::TraceLogger trace_log {};
        debug::print_genotype_log_likelihoods(stream(trace_log), genotypes,
                                              genotype_log_likilhoods, -1);
    } else if (DEBUG_MODE) {
        debug::print_genotype_log_likelihoods(stream(log), genotypes,
                                              genotype_log_likilhoods);
    }
    
    auto haplotype_prior_counts = compute_haplotype_prior_counts(haplotype_priors);
    
    //debug::print_haplotype_priors(haplotype_prior_counts, -1);
    
    const ModelConstants constants {haplotypes, genotypes, genotype_log_likilhoods, haplotype_prior_counts};
    
    auto haplotype_frequencies = init_haplotype_frequencies(haplotype_prior_counts, constants);
    
    haplotype_prior_counts.clear();
    
    auto genotype_log_marginals = init_genotype_log_marginals(genotypes, haplotype_frequencies);
    auto genotype_posteriors    = init_genotype_posteriors(genotype_log_marginals,
                                                           genotype_log_likilhoods);
    
    run_em(genotype_posteriors, haplotype_frequencies, genotype_log_marginals,
           constants, max_em_iterations_, em_epsilon_);
    
    return make_latents(haplotypes, std::move(genotypes), std::move(genotype_posteriors),
                        std::move(haplotype_frequencies), constants);
}

namespace debug
{
    void print_genotypes(const std::vector<Genotype<Haplotype>>& genotypes)
    {
        for (const auto& genotype : genotypes) {
            ::debug::print_variant_alleles(genotype);
        }
    }
    
    template <typename T, typename N>
    struct IsBigger
    {
        bool operator()(const std::pair<T, N>& lhs, const std::pair<T, N>& rhs) {
            return lhs.second > rhs.second;
        }
    };
    
    void print_haplotype_priors(const HaplotypePriorCountMap& prior_counts, const std::size_t n)
    {
        auto m = std::min(prior_counts.size(), n);
        
        std::cout << "Printing top " << m << " haplotype prior counts" << std::endl;
        
        std::vector<std::pair<Haplotype, double>> v {};
        v.reserve(prior_counts.size());
        
        std::copy(std::cbegin(prior_counts), std::cend(prior_counts), std::back_inserter(v));
        
        std::sort(std::begin(v), std::end(v), IsBigger<Haplotype, double>());
        
        for (unsigned i {0}; i < m; ++i) {
            ::debug::print_variant_alleles(v[i].first);
            std::cout << " " << std::setprecision(10) << v[i].second << std::endl;
        }
    }
    
    void print_genotype_log_marginals(const std::vector<Genotype<Haplotype>>& genotypes,
                                      const GenotypeLogMarginals& genotype_log_marginals,
                                      const std::size_t n)
    {
        print_genotype_log_marginals(std::cout, genotypes, genotype_log_marginals, n);
    }
    
    template <typename S>
    void print_genotype_log_likelihoods(S&&stream,
                                        const std::vector<Genotype<Haplotype>>& genotypes,
                                        const GenotypeLogLikelihoodMap& log_likelihoods,
                                        std::size_t n)
    {
        const auto m = std::min(genotypes.size(), n);
        
        if (m == genotypes.size()) {
            stream << "Printing all genotype likelihoods for each sample" << '\n';
        } else {
            stream << "Printing top " << m << " genotype likelihoods for each sample" << '\n';
        }
        
        for (const auto& sample_likelihoods : log_likelihoods) {
            stream << sample_likelihoods.first << ":" << '\n';
            
            std::vector<std::pair<Genotype<Haplotype>, double>> likelihoods {};
            likelihoods.reserve(sample_likelihoods.second.size());
            
            std::transform(std::cbegin(genotypes), std::cend(genotypes),
                           std::cbegin(sample_likelihoods.second),
                           std::back_inserter(likelihoods),
                           [] (const auto& genotype, auto log_liklihood) {
                               return std::make_pair(genotype, log_liklihood);
                           });
            
            const auto mth = std::next(std::begin(likelihoods), m);
            
            std::partial_sort(std::begin(likelihoods), mth, std::end(likelihoods),
                              [] (const auto& lhs, const auto& rhs) {
                                  return lhs.second > rhs.second;
                              });
            
            std::for_each(std::begin(likelihoods), mth,
                          [&] (const auto& p) {
                              ::debug::print_variant_alleles(stream, p.first);
                              stream << " " << std::setprecision(10) << p.second << '\n';
                          });
        }
    }
    
    void print_genotype_log_likelihoods(const std::vector<Genotype<Haplotype>>& genotypes,
                                        const GenotypeLogLikelihoodMap& log_likelihoods,
                                        const std::size_t n)
    {
        print_genotype_log_likelihoods(std::cout, genotypes, log_likelihoods, n);
    }
    
    template <typename S>
    void print_haplotype_frequencies(S&& stream,
                                     const HaplotypeFrequencyMap& haplotype_frequencies,
                                     const std::size_t n)
    {
        const auto m = std::min(haplotype_frequencies.size(), n);
        
        if (m == haplotype_frequencies.size()) {
            stream << "Printing all haplotype frequencies" << '\n';
        } else {
            stream << "Printing top " << m << " haplotype frequencies" << '\n';
        }
        
        std::vector<std::pair<Haplotype, double>> v {};
        v.reserve(haplotype_frequencies.size());
        
        std::copy(std::cbegin(haplotype_frequencies), std::cend(haplotype_frequencies),
                  std::back_inserter(v));
        
        const auto mth = std::next(std::begin(v), m);
        
        std::partial_sort(std::begin(v), mth, std::end(v),
                          [] (const auto& lhs, const auto& rhs) {
                              return lhs.second > rhs.second;
                          });
        
        std::for_each(std::begin(v), mth,
                      [&stream] (const auto& p) {
                          ::debug::print_variant_alleles(stream, p.first);
                          stream << " " << p.second << '\n';
                      });
    }
    
    void print_haplotype_frequencies(const HaplotypeFrequencyMap& haplotype_frequencies,
                                     const std::size_t n)
    {
        print_haplotype_frequencies(std::cout, haplotype_frequencies, n);
    }
    
    template <typename S>
    void print_genotype_log_marginals(S&& stream,
                                      const std::vector<Genotype<Haplotype>>& genotypes,
                                      const GenotypeLogMarginals& genotype_log_marginals,
                                      const std::size_t n)
    {
        const auto m = std::min(genotypes.size(), n);
        
        if (m == genotypes.size()) {
            stream << "Printing all genotype log marginals" << '\n';
        } else {
            stream << "Printing top " << m << " genotype log marginals" << '\n';
        }
        
        std::vector<std::pair<Genotype<Haplotype>, double>> v {};
        v.reserve(genotypes.size());
        
        std::transform(std::cbegin(genotypes), std::cend(genotypes),
                       std::cbegin(genotype_log_marginals),
                       std::back_inserter(v),
                       [] (const auto& genotype, auto log_marginal) {
                           return std::make_pair(genotype, log_marginal.log_probability);
                       });
        
        const auto mth = std::next(std::begin(v), m);
        
        std::partial_sort(std::begin(v), mth, std::end(v),
                          [] (const auto& lhs, const auto& rhs) {
                              return lhs.second > rhs.second;
                          });
        
        std::for_each(std::begin(v), mth,
                      [&stream] (const auto& p) {
                          ::debug::print_variant_alleles(stream, p.first);
                          stream << " " << p.second << '\n';
                      });
    }
    
    template <typename S>
    void print_genotype_posteriors(S&& stream,
                                   const std::vector<Genotype<Haplotype>>& genotypes,
                                   const GenotypePosteriorMap& genotype_posteriors,
                                   const std::size_t n)
    {
        const auto m = std::min(genotypes.size(), n);
        
        if (m == genotypes.size()) {
            stream << "Printing all genotype posterior for each sample" << '\n';
        } else {
            stream << "Printing top " << m << " genotype posterior for each sample" << '\n';
        }
        
        for (const auto& sample_posteriors : genotype_posteriors) {
            stream << sample_posteriors.first << ":" << '\n';
            
            std::vector<std::pair<Genotype<Haplotype>, double>> posteriors {};
            posteriors.reserve(sample_posteriors.second.size());
            
            std::transform(std::cbegin(genotypes), std::cend(genotypes),
                           std::cbegin(sample_posteriors.second),
                           std::back_inserter(posteriors),
                           [] (const auto& genotype, auto posterior) {
                               return std::make_pair(genotype, posterior);
                           });
            
            const auto mth = std::next(std::begin(posteriors), m);
            
            std::partial_sort(std::begin(posteriors), mth, std::end(posteriors),
                              [] (const auto& lhs, const auto& rhs) {
                                  return lhs.second > rhs.second;
                              });
            
            std::for_each(std::begin(posteriors), mth,
                          [&stream] (const auto& p) {
                              ::debug::print_variant_alleles(stream, p.first);
                              stream << " " << p.second << '\n';
                          });
        }
    }
    
    void print_genotype_posteriors(const std::vector<Genotype<Haplotype>>& genotypes,
                                   const GenotypePosteriorMap& genotype_posteriors,
                                   const std::size_t n)
    {
        print_genotype_posteriors(std::cout, genotypes, genotype_posteriors, n);
    }
} // namespace debug
} // namespace GenotypeModel
} // namespace Octopus

