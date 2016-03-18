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
    
    struct EMConstants
    {
        using InverseGenotypeTable = std::vector<std::vector<std::size_t>>;
        
        double norm;
        InverseGenotypeTable genotypes_containing_haplotypes;
        
        EMConstants(const double norm, InverseGenotypeTable&& table)
        : norm {norm}, genotypes_containing_haplotypes {std::move(table)} {}
    };
    } // namespace
    
    namespace debug
    {
        void print_genotypes(const std::vector<Genotype<Haplotype>>& genotypes);
        void print_haplotype_priors(const HaplotypePriorCountMap& prior_counts,
                                    std::size_t n = 5);
        void print_haplotype_frequencies(const HaplotypeFrequencyMap& haplotype_frequencies,
                                         std::size_t n = 5);
        template <typename S>
        void print_genotype_log_likelihoods(S&&stream,
                                            const std::vector<Genotype<Haplotype>>& genotypes,
                                            const GenotypeLogLikelihoodMap& log_likelihoods,
                                            std::size_t n = 5);
        void print_genotype_log_likelihoods(const std::vector<Genotype<Haplotype>>& genotypes,
                                            const GenotypeLogLikelihoodMap& log_likelihoods,
                                            std::size_t n = 5);
        void print_genotype_log_marginals(const std::vector<Genotype<Haplotype>>& genotypes,
                                          const GenotypeLogMarginals& genotype_log_marginals,
                                          std::size_t n = 10);
        void print_genotype_posteriors(const std::vector<Genotype<Haplotype>>& genotypes,
                                       const GenotypePosteriorMap& genotype_posteriors,
                                       std::size_t n = 10);
    } // namespace debug
    
    // non member methods
    
    namespace
    {
    HaplotypeFrequencyMap
    init_haplotype_frequencies(const HaplotypePriorCountMap& haplotype_prior_counts,
                               const double prior_count_sum)
    {
        HaplotypeFrequencyMap result {haplotype_prior_counts.size()};
        
        for (const auto& haplotype_count : haplotype_prior_counts) {
            result.emplace(haplotype_count.first, haplotype_count.second / prior_count_sum);
        }
        
        return result;
    }
    
    auto flatten_haplotype_prior_counts(const std::vector<Haplotype>& haplotypes,
                                        const HaplotypePriorCountMap& prior_counts)
    {
        std::vector<double> result(haplotypes.size());
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
                                         std::forward_as_tuple(std::ref(haplotype)),
                                         std::forward_as_tuple());
            it.first->second.reserve(r);
        }
        
        for (std::size_t i {0}; i < genotypes.size(); ++i) {
            for (const auto& haplotype : genotypes[i]) {
                result_map.at(haplotype).emplace_back(i);
            }
        }
        
        EMConstants::InverseGenotypeTable result {};
        
        for (const auto& haplotype : haplotypes) {
            result.emplace_back(std::move(result_map.at(haplotype)));
        }
        
        return result;
    }
    
    EMConstants
    make_haplotype_em_constants(const std::vector<Haplotype>& haplotypes,
                                const std::vector<Genotype<Haplotype>>& genotypes,
                                const double& prior_count_sum)
    {
        const auto norm = genotypes.size() * genotypes.front().ploidy() + prior_count_sum;
        return EMConstants {norm, make_inverse_genotype_table(haplotypes, genotypes)};
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
        using std::cbegin; using std::cend; using std::begin; using std::transform;
        
        GenotypePosteriorMap result {genotype_log_likilhoods.size()};
        
        for (const auto& sample_genotype_log_likilhoods : genotype_log_likilhoods) {
            SampleGenotypePosteriors sample_result(genotype_log_marginals.size());
            
            transform(cbegin(genotype_log_marginals), cend(genotype_log_marginals),
                      cbegin(sample_genotype_log_likilhoods.second),
                      begin(sample_result),
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
        using std::cbegin; using std::cend; using std::begin; using std::transform;
        
        for (auto& sample_genotype_posteriors : current_genotype_posteriors) {
            transform(cbegin(genotype_log_marginals), cend(genotype_log_marginals),
                      cbegin(genotype_log_likilhoods.at(sample_genotype_posteriors.first)),
                      begin(sample_genotype_posteriors.second),
                      [] (const auto& genotype_log_marginal, const auto& genotype_log_likilhood) {
                          return genotype_log_marginal.log_probability + genotype_log_likilhood;
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
                           std::cbegin(sample_posteriors.second), std::begin(result),
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
                                        const EMConstants& constants)
    {
        const auto collaped_posteriors = collapse_genotype_posteriors(genotype_posteriors);
        
        double max_frequency_change {0};
        
        for (std::size_t i {0}; i < haplotypes.size(); ++i) {
            auto& current_frequency = current_haplotype_frequencies.at(haplotypes[i]);
            
            double new_frequency {haplotype_prior_counts[i]};
            
            for (const auto& genotype_index : constants.genotypes_containing_haplotypes[i]) {
                new_frequency += collaped_posteriors[genotype_index];
            }
            
            new_frequency /= constants.norm;
            
            const auto frequency_change = std::abs(current_frequency - new_frequency);
            
            if (frequency_change > max_frequency_change) {
                max_frequency_change = frequency_change;
            }
            
            current_frequency = new_frequency;
        }
        
        return max_frequency_change;
    }
    
    double do_em_iteration(const std::vector<Haplotype>& haplotypes,
                           GenotypePosteriorMap& genotype_posteriors,
                           HaplotypeFrequencyMap& haplotype_frequencies,
                           GenotypeLogMarginals& genotype_log_marginals,
                           const GenotypeLogLikelihoodMap& genotype_log_likilhoods,
                           const HaplotypePriorCounts& haplotype_prior_counts,
                           const EMConstants& em_constants)
    {
        const auto max_change = update_haplotype_frequencies(haplotypes,
                                                             haplotype_frequencies,
                                                             haplotype_prior_counts,
                                                             genotype_posteriors,
                                                             em_constants);
        
        update_genotype_log_marginals(genotype_log_marginals, haplotype_frequencies);
        
        update_genotype_posteriors(genotype_posteriors, genotype_log_marginals, genotype_log_likilhoods);
        
        return max_change;
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
        
        Population::Latents::HaplotypeFrequencyMap haps {{genotype[0], 1}};
        
        return Population::Latents {std::move(result), std::move(haps)};
    }
    
    Population::Latents
    make_latents(std::vector<Genotype<Haplotype>>&& genotypes,
                 GenotypePosteriorMap&& genotype_posteriors,
                 HaplotypeFrequencyMap&& haplotype_frequencies)
    {
        Population::Latents::GenotypeProbabilityMap result {std::make_move_iterator(std::begin(genotypes)),
                                                    std::make_move_iterator(std::end(genotypes))};
        
        insert_samples(std::move(genotype_posteriors), result);
        
        return Population::Latents {std::move(result), std::move(haplotype_frequencies)};
    }
    } // namespace
    
    // private methods
    
    Population::Latents
    Population::infer_latents(const std::vector<SampleIdType>& samples,
                              const std::vector<Haplotype>& haplotypes,
                              const HaplotypePrioMap& haplotype_priors,
                              const HaplotypeLikelihoodCache& haplotype_likelihoods)
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
        
        if (DEBUG_MODE) {
            debug::print_genotype_log_likelihoods(stream(log), genotypes, genotype_log_likilhoods, -1);
        }
        
        auto haplotype_prior_count_map = compute_haplotype_prior_counts(haplotype_priors);
        
        auto haplotype_prior_counts    = flatten_haplotype_prior_counts(haplotypes,
                                                                        haplotype_prior_count_map);
        const auto prior_count_sum     = std::accumulate(std::cbegin(haplotype_prior_counts),
                                                         std::cend(haplotype_prior_counts), 0.0);
        
        const auto em_constants = make_haplotype_em_constants(haplotypes, genotypes, prior_count_sum);
        
        auto haplotype_frequencies = init_haplotype_frequencies(haplotype_prior_count_map, prior_count_sum);
        
        haplotype_prior_count_map.clear();
        
        auto genotype_log_marginals = init_genotype_log_marginals(genotypes, haplotype_frequencies);
        auto genotype_posteriors    = init_genotype_posteriors(genotype_log_marginals, genotype_log_likilhoods);
        
        for (unsigned n {0}; n < max_em_iterations_; ++n) {
            if (do_em_iteration(haplotypes, genotype_posteriors, haplotype_frequencies,
                                genotype_log_marginals, genotype_log_likilhoods,
                                haplotype_prior_counts, em_constants) < em_epsilon_) break;
        }
        
        return make_latents(std::move(genotypes), std::move(genotype_posteriors),
                            std::move(haplotype_frequencies));
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
            
            std::cout << "printing top " << m << " haplotype prior counts" << std::endl;
            
            std::vector<std::pair<Haplotype, double>> v {};
            v.reserve(prior_counts.size());
            
            std::copy(std::cbegin(prior_counts), std::cend(prior_counts), std::back_inserter(v));
            
            std::sort(std::begin(v), std::end(v), IsBigger<Haplotype, double>());
            
            for (unsigned i {0}; i < m; ++i) {
                ::debug::print_variant_alleles(v[i].first);
                std::cout << " " << std::setprecision(10) << v[i].second << std::endl;
            }
        }
        
        void print_haplotype_frequencies(const HaplotypeFrequencyMap& haplotype_frequencies,
                                         const std::size_t n)
        {
            auto m = std::min(haplotype_frequencies.size(), n);
            
            std::cout << "printing top " << m << " haplotype frequencies" << std::endl;
            
            std::vector<std::pair<Haplotype, double>> v {};
            v.reserve(haplotype_frequencies.size());
            
            std::copy(std::cbegin(haplotype_frequencies), std::cend(haplotype_frequencies), std::back_inserter(v));
            
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
            auto m = std::min(genotypes.size(), n);
            
            std::cout << "printing top " << m << " genotype log marginals" << std::endl;
            
            std::vector<std::pair<Genotype<Haplotype>, double>> v {};
            v.reserve(genotypes.size());
            
            std::transform(std::cbegin(genotypes), std::cend(genotypes),
                           std::cbegin(genotype_log_marginals),
                           std::back_inserter(v),
                           [] (const auto& genotype, auto log_marginal) {
                               return std::make_pair(genotype, log_marginal.log_probability);
                           });
            
            std::sort(std::begin(v), std::end(v), IsBigger<Genotype<Haplotype>, double>());
            
            for (unsigned i {0}; i < std::min(n, v.size()); ++i) {
                ::debug::print_variant_alleles(v[i].first);
                std::cout << " " << std::setprecision(10) << v[i].second << std::endl;
            }
        }
        
        template <typename S>
        void print_genotype_log_likelihoods(S&&stream,
                                            const std::vector<Genotype<Haplotype>>& genotypes,
                                            const GenotypeLogLikelihoodMap& log_likelihoods,
                                            std::size_t n)
        {
            const auto m = std::min(genotypes.size(), n);
            
            if (m == genotypes.size()) {
                stream << "printing all genotype likelihoods for each sample" << '\n';
            } else {
                stream << "printing top " << m << " genotype likelihoods for each sample" << '\n';
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
        
        void print_genotype_posteriors(const std::vector<Genotype<Haplotype>>& genotypes,
                                       const GenotypePosteriorMap& genotype_posteriors,
                                       const std::size_t n)
        {
            const auto m = std::min(genotypes.size(), n);
            
            std::cout << "printing top " << m << " genotype posterior for each sample" << std::endl;
            
            for (const auto& sample_posteriors : genotype_posteriors) {
                std::cout << sample_posteriors.first << ":" << std::endl;
                
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
                              [] (const auto& p) {
                                  ::debug::print_variant_alleles(p.first);
                                  std::cout << " " << std::setprecision(10) << p.second << '\n';
                              });
            }
        }
    } // namespace debug
    } // namespace GenotypeModel
} // namespace Octopus

