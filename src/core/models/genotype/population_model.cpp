//
//  population_model.cpp
//  Octopus
//
//  Created by Daniel Cooke on 26/08/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "population_model.hpp"

#include <utility>
#include <algorithm>
#include <cmath>
#include <cassert>
#include <iostream>

#include <utils/maths.hpp>

#include "germline_likelihood_model.hpp"

namespace octopus { namespace model {

PopulationModel::PopulationModel(const CoalescentModel& genotype_prior_model)
:
genotype_prior_model_ {genotype_prior_model}
{}

namespace
{

using GenotypeLogLikelihoodVector  = std::vector<double>;
using SampleGenotypeLogLikelihoods = std::vector<GenotypeLogLikelihoodVector>;

struct GenotypeLogProbability
{
    const Genotype<Haplotype>& genotype;
    double log_probability;
};

using GenotypeLogMarginals = std::vector<GenotypeLogProbability>;

using GenotypePosteriorVector  = std::vector<double>;
using SampleGenotypePosteriors = std::vector<GenotypePosteriorVector>;

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

using HaplotypeReference    = std::reference_wrapper<const Haplotype>;
using HaplotypeFrequencyMap = std::unordered_map<HaplotypeReference, double>;

double calculate_frequency_update_norm(const std::size_t num_samples, const unsigned ploidy)
{
    return static_cast<double>(num_samples) * ploidy;
}

struct ModelConstants
{
    using InverseGenotypeTable = std::vector<std::vector<std::size_t>>;
    
    const std::vector<Haplotype>& haplotypes;
    const std::vector<Genotype<Haplotype>>& genotypes;
    const SampleGenotypeLogLikelihoods& genotype_log_likilhoods;
    const unsigned ploidy;
    const double frequency_update_norm;
    const InverseGenotypeTable genotypes_containing_haplotypes;
};

} // namespace

namespace debug {
//        void print_genotypes(const std::vector<Genotype<Haplotype>>& genotypes);
//        template <typename S>
//        void print_genotype_log_likelihoods(S&&stream,
//                                            const std::vector<Genotype<Haplotype>>& genotypes,
//                                            const SampleGenotypeLogLikelihoods& log_likelihoods,
//                                            std::size_t n = 5);
//        void print_genotype_log_likelihoods(const std::vector<Genotype<Haplotype>>& genotypes,
//                                            const SampleGenotypeLogLikelihoods& log_likelihoods,
//                                            std::size_t n = 5);
//        template <typename S>
//        void print_haplotype_frequencies(S&& stream,
//                                         const HaplotypeFrequencyMap& haplotype_frequencies,
//                                         const std::size_t n = 5);
//        void print_haplotype_frequencies(const HaplotypeFrequencyMap& haplotype_frequencies,
//                                         std::size_t n = 5);
//        template <typename S>
//        void print_genotype_log_marginals(S&& stream,
//                                          const std::vector<Genotype<Haplotype>>& genotypes,
//                                          const GenotypeLogMarginals& genotype_log_marginals,
//                                          std::size_t n = 5);
//        void print_genotype_log_marginals(const std::vector<Genotype<Haplotype>>& genotypes,
//                                          const GenotypeLogMarginals& genotype_log_marginals,
//                                          std::size_t n = 5);
//        template <typename S>
//        void print_genotype_posteriors(S&& stream,
//                                       const std::vector<Genotype<Haplotype>>& genotypes,
//                                       const GenotypePosteriorMap& genotype_posteriors,
//                                       std::size_t n = 5);
//        void print_genotype_posteriors(const std::vector<Genotype<Haplotype>>& genotypes,
//                                       const GenotypePosteriorMap& genotype_posteriors,
//                                       std::size_t n = 5);
} // namespace debug

// non member methods

namespace
{

HaplotypeFrequencyMap
init_haplotype_frequencies(const ModelConstants& constants)
{
    HaplotypeFrequencyMap result {constants.haplotypes.size()};
    
    for (const auto& haplotype : constants.haplotypes) {
        result.emplace(haplotype, 1.0 / constants.haplotypes.size());
    }
    
    return result;
}

namespace detail
{
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
        
        static const double ln_2 {std::log(2.0)};
        
        return std::log(haplotype_frequencies.at(genotype[0]))
        + std::log(haplotype_frequencies.at(genotype[1])) + ln_2;
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
}

// TODO: improve this, possible bottleneck in EM update at the moment
template <typename Genotype, typename Map>
double log_hardy_weinberg(const Genotype& genotype, const Map& haplotype_frequencies)
{
    switch (genotype.ploidy()) {
        case 1 : return detail::log_hardy_weinberg_haploid(genotype, haplotype_frequencies);
        case 2 : return detail::log_hardy_weinberg_diploid(genotype, haplotype_frequencies);
        case 3 : return detail::log_hardy_weinberg_triploid(genotype, haplotype_frequencies);
        default: return detail::log_hardy_weinberg_polyploid(genotype, haplotype_frequencies);
    }
}

SampleGenotypeLogLikelihoods
compute_genotype_log_likelihoods(const std::vector<SampleName>& samples,
                                 const std::vector<Genotype<Haplotype>>& genotypes,
                                 const HaplotypeLikelihoodCache& haplotype_likelihoods)
{
    assert(!genotypes.empty());
    
    GermlineLikelihoodModel likelihood_model {haplotype_likelihoods};
    
    SampleGenotypeLogLikelihoods result {};
    result.reserve(samples.size());
    
    std::transform(std::cbegin(samples), std::cend(samples), std::back_inserter(result),
                   [&genotypes, &haplotype_likelihoods, &likelihood_model] (const auto& sample) {
                       GenotypeLogLikelihoodVector likelihoods(genotypes.size());
                       
                       haplotype_likelihoods.prime(sample);
                       
                       std::transform(std::cbegin(genotypes), std::cend(genotypes),
                                      std::begin(likelihoods),
                                      [&sample, &likelihood_model] (const auto& genotype) {
                                          return likelihood_model.ln_likelihood(genotype);
                                      });
                       
                       return likelihoods;
                   });
    
    return result;
}

GenotypeLogMarginals
init_genotype_log_marginals(const std::vector<Genotype<Haplotype>>& genotypes,
                            const HaplotypeFrequencyMap& haplotype_frequencies)
{
    GenotypeLogMarginals result {};
    result.reserve(genotypes.size());
    
    for (const auto& genotype : genotypes) {
        result.push_back({genotype, log_hardy_weinberg(genotype, haplotype_frequencies)});
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

SampleGenotypePosteriors
init_genotype_posteriors(const GenotypeLogMarginals& genotype_log_marginals,
                         const SampleGenotypeLogLikelihoods& genotype_log_likilhoods)
{
    SampleGenotypePosteriors result {};
    result.reserve(genotype_log_likilhoods.size());
    
    for (const auto& sample_genotype_log_likilhoods : genotype_log_likilhoods) {
        GenotypePosteriorVector posteriors(genotype_log_marginals.size());
        
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

void update_genotype_posteriors(SampleGenotypePosteriors& current_genotype_posteriors,
                                const GenotypeLogMarginals& genotype_log_marginals,
                                const SampleGenotypeLogLikelihoods& genotype_log_likilhoods)
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

auto collapse_genotype_posteriors(const SampleGenotypePosteriors& genotype_posteriors)
{
    assert(!genotype_posteriors.empty());
    
    std::vector<double> result(genotype_posteriors.front().size());
    
    for (const auto& sample_posteriors : genotype_posteriors) {
        std::transform(std::cbegin(result), std::cend(result),
                       std::cbegin(sample_posteriors),
                       std::begin(result),
                       [] (const auto curr, const auto p) {
                           return curr + p;
                       });
    }
    
    return result;
}

double update_haplotype_frequencies(const std::vector<Haplotype>& haplotypes,
                                    HaplotypeFrequencyMap& current_haplotype_frequencies,
                                    const SampleGenotypePosteriors& genotype_posteriors,
                                    const ModelConstants::InverseGenotypeTable& genotypes_containing_haplotypes,
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

double do_em_iteration(SampleGenotypePosteriors& genotype_posteriors,
                       HaplotypeFrequencyMap& haplotype_frequencies,
                       GenotypeLogMarginals& genotype_log_marginals,
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

void run_em(SampleGenotypePosteriors& genotype_posteriors,
            HaplotypeFrequencyMap& haplotype_frequencies,
            GenotypeLogMarginals& genotype_log_marginals,
            const ModelConstants& constants,
            const unsigned max_iterations, const double epsilon)
{
//            if (TRACE_MODE) {
//                logging::TraceLogger log {};
//                stream(log) << "Running EM";
//                debug::print_haplotype_frequencies(stream(log), haplotype_frequencies, -1);
//                debug::print_genotype_log_marginals(stream(log), constants.genotypes, genotype_log_marginals, -1);
//                //        debug::print_genotype_posteriors(stream(log), constants.genotypes, genotype_posteriors, -1);
//            }
    
    for (unsigned n {1}; n <= max_iterations; ++n) {
        const auto max_change = do_em_iteration(genotype_posteriors, haplotype_frequencies,
                                                genotype_log_marginals,constants);
        
//                if (TRACE_MODE) {
//                    logging::TraceLogger log {};
//                    stream(log) << "EM iteration " << n << " : max change = " << max_change;
//                    debug::print_haplotype_frequencies(stream(log), haplotype_frequencies, -1);
//                    //            debug::print_genotype_log_marginals(stream(log), constants.genotypes, genotype_log_marginals, -1);
//                    //            debug::print_genotype_posteriors(stream(log), constants.genotypes, genotype_posteriors, -1);
//                }
        
        if (max_change <= epsilon) break;
    }
}

//        Population::Latents
//        make_single_genotype_latents(const std::vector<SampleName>& samples,
//                                     const Genotype<Haplotype>& genotype)
//        {
//            Population::Latents::GenotypeProbabilityMap result {};
//            result.push_back(genotype);
//            
//            for (const auto& sample : samples) {
//                insert_sample(sample, std::vector<double> {1}, result);
//            }
//            
//            Population::Latents::HaplotypeFrequencyMap haplotype_frequencies {{genotype[0], 1.0}};
//            Population::Latents::HaplotypeProbabilityMap haplotype_posteriors {{genotype[0], 1.0}};
//            
//            return Population::Latents {
//                std::move(haplotype_frequencies),
//                std::move(result),
//                std::move(haplotype_posteriors)
//            };
//        }

auto calculate_haplotype_posteriors(const std::vector<Haplotype>& haplotypes,
                                    const std::vector<Genotype<Haplotype>>& genotypes,
                                    const GenotypePosteriorVector& genotype_posteriors)
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
                                    const SampleGenotypePosteriors& genotype_posteriors,
                                    const ModelConstants::InverseGenotypeTable& inverse_genotypes)
{
    if (genotype_posteriors.size() == 1) {
        // can optimise this case
        return calculate_haplotype_posteriors(haplotypes, genotypes, genotype_posteriors.front());
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
                                                     return curr + sample_genotype_posteriors[i];
                                                 });
        }
        
        result.emplace(haplotype, 1.0 - prob_not_observed);
        
        ++it;
    }
    
    return result;
}

PopulationModel::InferredLatents
make_latents(const std::vector<Haplotype>& haplotypes,
             const std::vector<Genotype<Haplotype>>& genotypes,
             SampleGenotypePosteriors&& genotype_posteriors,
             const ModelConstants& constants)
{
    auto haplotype_posteriors = calculate_haplotype_posteriors(haplotypes, genotypes,
                                                               genotype_posteriors,
                                                               constants.genotypes_containing_haplotypes);
    
    return PopulationModel::InferredLatents {
        std::move(genotype_posteriors),
        std::move(haplotype_posteriors),
        1.0 // TODO
    };
}

} // namespace

PopulationModel::InferredLatents
PopulationModel::infer_latents(const std::vector<SampleName>& samples, const GenotypeVector& genotypes,
                               const std::vector<Haplotype>& haplotypes,
                               const HaplotypeLikelihoodCache& haplotype_likelihoods) const
{
    assert(!genotypes.empty());
    
//    if (genotypes.size() == 1) {
//        return make_single_genotype_latents(samples, genotypes.front());
//    }
    
    const auto genotype_log_likilhoods = compute_genotype_log_likelihoods(samples, genotypes,
                                                                          haplotype_likelihoods);
    
//    if (TRACE_MODE) {
//        logging::TraceLogger trace_log {};
//        debug::print_genotype_log_likelihoods(stream(trace_log), genotypes,
//                                              genotype_log_likilhoods, -1);
//    } else if (DEBUG_MODE) {
//        logging::DebugLogger log {};
//        debug::print_genotype_log_likelihoods(stream(log), genotypes,
//                                              genotype_log_likilhoods);
//    }
    
    const ModelConstants constants {haplotypes, genotypes, genotype_log_likilhoods};
    
    auto haplotype_frequencies = init_haplotype_frequencies(constants);
    
    auto genotype_log_marginals = init_genotype_log_marginals(genotypes, haplotype_frequencies);
    auto genotype_posteriors    = init_genotype_posteriors(genotype_log_marginals,
                                                           genotype_log_likilhoods);
    
    run_em(genotype_posteriors, haplotype_frequencies, genotype_log_marginals,
           constants, 100, 0.001);
    
    return make_latents(haplotypes, genotypes, std::move(genotype_posteriors), constants);
}

namespace debug {
    
} // namespace debug
} // namespace model
} // namespace octopus

