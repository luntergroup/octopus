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

#include "hardy_weinberg_model.hpp"
#include "dirichlet_model.hpp"
#include "fixed_ploidy_genotype_likelihood_model.hpp"
#include "maths.hpp"
#include "read_utils.hpp"

#include <iostream> // TEST
#include <chrono>   // TEST

namespace Octopus
{
    namespace GenotypeModel
    {
    // public methods
    
    Population::Population(unsigned ploidy, unsigned max_em_iterations, double em_epsilon)
    :
    ploidy_ {ploidy},
    max_em_iterations_ {max_em_iterations},
    em_epsilon_ {em_epsilon}
    {}
    
    // non member methods
    
    namespace {
    
    using GenotypeLogLikelihood         = double;
    using SampleGenotypeLogLikelihoods  = std::vector<GenotypeLogLikelihood>;
    using GenotypeLogLikelihoodMap      = std::unordered_map<SampleIdType, SampleGenotypeLogLikelihoods>;
    
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
    using GenotypePosteriors       = std::unordered_map<SampleIdType, SampleGenotypePosteriors>;
    
    using HaplotypeLogFrequencies = std::unordered_map<Haplotype, double>;
    } // namespace
    
    namespace debug {
        void print_haplotype_priors(const HaplotypePriorCountMap& prior_counts, size_t n = 5);
        void print_haplotype_frequencies(const HaplotypeFrequencyMap& haplotype_frequencies, size_t n = 5);
        void print_genotype_log_likelihoods(const std::vector<Genotype<Haplotype>>& genotypes,
                                            const GenotypeLogLikelihoodMap& log_likelihoods, size_t n = 5);
        void print_genotype_log_marginals(const std::vector<Genotype<Haplotype>>& genotypes,
                                          const GenotypeLogMarginals& genotype_log_marginals,
                                          size_t n = 10);
        void print_genotype_posteriors(const std::vector<Genotype<Haplotype>>& genotypes,
                                       const GenotypePosteriors& genotype_posteriors,
                                       size_t n = 10);
    } // namespace debug
    
    // non member methods
    
    namespace
    {
    GenotypeLogLikelihoodMap
    compute_genotype_log_likelihoods(const std::vector<Genotype<Haplotype>>& genotypes,
                                     const ReadMap& reads,
                                     HaplotypeLikelihoodCache& haplotype_likelihoods)
    {
        GenotypeLogLikelihoodMap result {};
        
        if (genotypes.empty()) return result;
        
        FixedPloidyGenotypeLikelihoodModel read_model {genotypes.front().ploidy(), haplotype_likelihoods};
        
        result.reserve(reads.size());
        
        for (const auto& sample_reads_p : reads) {
            const auto& sample       = sample_reads_p.first;
            const auto& sample_reads = sample_reads_p.second;
            
            SampleGenotypeLogLikelihoods sample_log_likelihoods(genotypes.size());
            
            std::transform(std::cbegin(genotypes), cend(genotypes), std::begin(sample_log_likelihoods),
                           [&sample_reads, &read_model] (const auto& genotype) {
                               return log_probability(sample_reads, genotype, read_model);
                           });
            
            result.emplace(sample, std::move(sample_log_likelihoods));
        }
        
        haplotype_likelihoods.clear();
        
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
        const auto norm = Maths::log_sum_exp<double>(unnormalised_log_probabilities);
        
        std::transform(std::cbegin(unnormalised_log_probabilities),
                       std::cend(unnormalised_log_probabilities),
                       std::begin(unnormalised_log_probabilities),
                       [norm] (const double p) {
                           return std::exp(p - norm);
                       });
    }
    
    GenotypePosteriors
    init_genotype_posteriors(const GenotypeLogMarginals& genotype_log_marginals,
                             const GenotypeLogLikelihoodMap& genotype_log_likilhoods)
    {
        using std::cbegin; using std::cend; using std::begin; using std::transform;
        
        GenotypePosteriors result {};
        result.reserve(genotype_log_likilhoods.size());
        
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
    
    void update_genotype_posteriors(GenotypePosteriors& current_genotype_posteriors,
                                    const HaplotypeFrequencyMap& haplotype_frequencies,
                                    const GenotypeLogMarginals& genotype_log_marginals,
                                    const GenotypeLogLikelihoodMap& genotype_log_likilhoods)
    {
        using std::cbegin; using std::cend; using std::begin; using std::transform;
        
        for (auto& sample_genotype_posteriors : current_genotype_posteriors) {
            
            const auto& sample_genotype_log_likilhoods = genotype_log_likilhoods.at(sample_genotype_posteriors.first);
            
            transform(cbegin(genotype_log_marginals), cend(genotype_log_marginals),
                      cbegin(sample_genotype_log_likilhoods),
                      begin(sample_genotype_posteriors.second),
                      [] (const auto& genotype_log_marginal, const auto& genotype_log_likilhood) {
                          return genotype_log_marginal.log_probability + genotype_log_likilhood;
                      });
            
            normalise_exp(sample_genotype_posteriors.second);
        }
    }
    
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
    
    double update_haplotype_frequencies(HaplotypeFrequencyMap& current_haplotype_frequencies,
                                        const HaplotypePriorCountMap& haplotype_prior_counts,
                                        const std::vector<Genotype<Haplotype>>& genotypes,
                                        const GenotypePosteriors& genotype_posteriors,
                                        const double prior_count_sum)
    {
        using std::cbegin; using std::cend; using std::accumulate;
        
        double max_frequency_change {0};
        
        const double norm {genotype_posteriors.size() * genotypes.front().ploidy() + prior_count_sum};
        
        for (auto& haplotype_frequency : current_haplotype_frequencies) {
            const auto& haplotype = haplotype_frequency.first;
            
            double new_frequency {haplotype_prior_counts.at(haplotype)};
            
            for (size_t k {}; k < genotypes.size(); ++k) {
                const auto n = genotypes[k].count(haplotype);
                
                if (n > 0) {
                    auto s = accumulate(cbegin(genotype_posteriors), cend(genotype_posteriors), 0.0,
                                        [k] (double curr, const auto& sample_genotype_posteriors) {
                                            return curr + sample_genotype_posteriors.second[k];
                                        });
                    new_frequency += n * s;
                }
            }
            
            new_frequency /= norm;
            
            const auto frequency_change = std::abs(haplotype_frequency.second - new_frequency);
            
            if (frequency_change > max_frequency_change) {
                max_frequency_change = frequency_change;
            }
            
            haplotype_frequency.second = new_frequency;
        }
        
        return max_frequency_change;
    }
    
    double do_em_iteration(const std::vector<Haplotype>& haplotypes,
                           const std::vector<Genotype<Haplotype>>& genotypes,
                           GenotypePosteriors& genotype_posteriors,
                           HaplotypeFrequencyMap& haplotype_frequencies,
                           GenotypeLogMarginals& genotype_log_marginals,
                           const GenotypeLogLikelihoodMap& genotype_log_likilhoods,
                           const HaplotypePriorCountMap& haplotype_prior_counts,
                           const double prior_count_sum)
    {
        const auto max_change = update_haplotype_frequencies(haplotype_frequencies,
                                                             haplotype_prior_counts,
                                                             genotypes,
                                                             genotype_posteriors,
                                                             prior_count_sum);
        
        //debug::print_haplotype_frequencies(haplotype_frequencies);
        
        update_genotype_log_marginals(genotype_log_marginals, haplotype_frequencies);
        
        update_genotype_posteriors(genotype_posteriors, haplotype_frequencies,
                                   genotype_log_marginals, genotype_log_likilhoods);
        
        //debug::print_genotype_posteriors(genotypes, genotype_posteriors);
        
        return max_change;
    }
    
    Population::Latents
    make_single_genotype_latents(const Genotype<Haplotype>& genotype, const ReadMap& reads)
    {
        Population::Latents::GenotypeProbabilityMap result {};
        result.push_back(genotype);
        
        for (const auto& s : reads) {
            insert_sample(s.first, std::vector<double> {1}, result);
        }
        
        Population::Latents::HaplotypeFrequencyMap haps {{genotype.at(0), 1}};
        
        return Population::Latents {std::move(result), std::move(haps)};
    }
    
    Population::Latents
    make_latents(std::vector<Genotype<Haplotype>>&& genotypes,
                 GenotypePosteriors&& genotype_posteriors,
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
    Population::infer_latents(const std::vector<Haplotype>& haplotypes,
                              const HaplotypePrioMap& haplotype_priors,
                              HaplotypeLikelihoodCache& haplotype_likelihoods,
                              const ReadMap& reads)
    {
        assert(!haplotypes.empty());
        assert(!reads.empty());
        
        auto genotypes = generate_all_genotypes(haplotypes, ploidy_);
        
        //std::cout << "there are " << genotypes.size() << " candidate genotypes" << std::endl;
        
        if (genotypes.size() == 1) {
            return make_single_genotype_latents(genotypes.front(), reads);
        }
        
        const auto genotype_log_likilhoods = compute_genotype_log_likelihoods(genotypes, reads,
                                                                              haplotype_likelihoods);
        
        auto haplotype_prior_counts = compute_haplotype_prior_counts(haplotype_priors);
        const auto prior_count_sum  = Maths::sum_values(haplotype_prior_counts);
        
        auto haplotype_frequencies  = init_haplotype_frequencies(haplotype_prior_counts, prior_count_sum);
        auto genotype_log_marginals = init_genotype_log_marginals(genotypes, haplotype_frequencies);
        auto genotype_posteriors    = init_genotype_posteriors(genotype_log_marginals, genotype_log_likilhoods);
        
        for (unsigned n {}; n < max_em_iterations_; ++n) {
            //std::cout << "******* EM iteration " << n << " *******" << std::endl;
            if (do_em_iteration(haplotypes, genotypes, genotype_posteriors, haplotype_frequencies,
                                genotype_log_marginals, genotype_log_likilhoods,
                                haplotype_prior_counts, prior_count_sum) < em_epsilon_) break;
        }
        
        //std::cout << "finished EM" << std::endl;
        
        return make_latents(std::move(genotypes), std::move(genotype_posteriors), std::move(haplotype_frequencies));
    }
    
    namespace debug {
        template <typename T, typename N>
        struct IsBigger
        {
            bool operator()(const std::pair<T, N>& lhs, const std::pair<T, N>& rhs) {
                return lhs.second > rhs.second;
            }
        };
        
        void print_haplotype_priors(const HaplotypePriorCountMap& prior_counts, size_t n)
        {
            auto m = std::min(prior_counts.size(), n);
            
            std::cout << "printing top " << m << " haplotype prior counts" << std::endl;
            
            std::vector<std::pair<Haplotype, double>> v {};
            v.reserve(prior_counts.size());
            
            std::copy(std::cbegin(prior_counts), std::cend(prior_counts), std::back_inserter(v));
            
            std::sort(std::begin(v), std::end(v), IsBigger<Haplotype, double>());
            
            for (unsigned i {}; i < m; ++i) {
                print_variant_alleles(v[i].first);
                std::cout << " " << std::setprecision(10) << v[i].second << std::endl;
            }
        }
        
        void print_haplotype_frequencies(const HaplotypeFrequencyMap& haplotype_frequencies,
                                         size_t n)
        {
            auto m = std::min(haplotype_frequencies.size(), n);
            
            std::cout << "printing top " << m << " haplotype frequencies" << std::endl;
            
            std::vector<std::pair<Haplotype, double>> v {};
            v.reserve(haplotype_frequencies.size());
            
            std::copy(std::cbegin(haplotype_frequencies), std::cend(haplotype_frequencies), std::back_inserter(v));
            
            std::sort(std::begin(v), std::end(v), IsBigger<Haplotype, double>());
            
            for (unsigned i {}; i < m; ++i) {
                print_variant_alleles(v[i].first);
                std::cout << " " << std::setprecision(10) << v[i].second << std::endl;
            }
        }
        
        void print_genotype_log_marginals(const std::vector<Genotype<Haplotype>>& genotypes,
                                          const GenotypeLogMarginals& genotype_log_marginals,
                                          size_t n)
        {
            auto m = std::min(genotypes.size(), n);
            
            std::cout << "printing top " << m << " genotype log marginals" << std::endl;
            
            std::vector<std::pair<Genotype<Haplotype>, double>> v {};
            v.reserve(genotypes.size());
            
            std::transform(std::cbegin(genotypes), std::cend(genotypes),
                           std::cbegin(genotype_log_marginals), std::back_inserter(v),
                           [] (const auto& genotype, auto log_marginal) {
                               return std::make_pair(genotype, log_marginal.log_probability);
                           });
            
            std::sort(std::begin(v), std::end(v), IsBigger<Genotype<Haplotype>, double>());
            
            for (unsigned i {}; i < std::min(n, v.size()); ++i) {
                print_variant_alleles(v[i].first);
                std::cout << " " << std::setprecision(10) << v[i].second << std::endl;
            }
        }
        
        void print_genotype_log_likelihoods(const std::vector<Genotype<Haplotype>>& genotypes,
                                            const GenotypeLogLikelihoodMap& log_likelihoods,
                                            size_t n)
        {
            auto m = std::min(genotypes.size(), n);
            
            std::cout << "printing top " << m << " genotype likelihoods for each sample" << std::endl;
            
            for (const auto& sample_likelihoods : log_likelihoods) {
                std::cout << sample_likelihoods.first << ":" << std::endl;
                
                std::vector<std::pair<Genotype<Haplotype>, double>> v {};
                v.reserve(sample_likelihoods.second.size());
                
                std::transform(std::cbegin(genotypes), std::cend(genotypes),
                               std::cbegin(sample_likelihoods.second), std::back_inserter(v),
                               [] (const auto& genotype, auto log_liklihood) {
                                   return std::make_pair(genotype, log_liklihood);
                               });
                
                std::sort(std::begin(v), std::end(v), IsBigger<Genotype<Haplotype>, double>());
                
                for (unsigned i {}; i < std::min(n, v.size()); ++i) {
                    print_variant_alleles(v[i].first);
                    std::cout << " " << std::setprecision(10) << v[i].second << std::endl;
                }
            }
        }
        
        void print_genotype_posteriors(const std::vector<Genotype<Haplotype>>& genotypes,
                                       const GenotypePosteriors& genotype_posteriors,
                                       size_t n)
        {
            auto m = std::min(genotypes.size(), n);
            
            std::cout << "printing top " << m << " genotype posterior for each sample" << std::endl;
            
            for (const auto& sample_posteriors : genotype_posteriors) {
                std::cout << sample_posteriors.first << ":" << std::endl;
                
                std::vector<std::pair<Genotype<Haplotype>, double>> v {};
                v.reserve(sample_posteriors.second.size());
                
                std::transform(std::cbegin(genotypes), std::cend(genotypes),
                               std::cbegin(sample_posteriors.second), std::back_inserter(v),
                               [] (const auto& genotype, auto posterior) {
                                   return std::make_pair(genotype, posterior);
                               });
                
                std::sort(std::begin(v), std::end(v), IsBigger<Genotype<Haplotype>, double>());
                
                for (unsigned i {}; i < std::min(n, v.size()); ++i) {
                    print_variant_alleles(v[i].first);
                    std::cout << " " << std::setprecision(10) << v[i].second << std::endl;
                }
            }
        }
    } // namespace debug
    } // namespace GenotypeModel
} // namespace Octopus

