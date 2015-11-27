//
//  population_genotype_model.cpp
//  Octopus
//
//  Created by Daniel Cooke on 26/08/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "population_genotype_model.hpp"

#include <algorithm> // std::transform, std::for_each
#include <numeric>   // std::accumulate
#include <cmath>

#include "read_model.hpp"
#include "maths.hpp"
#include "haplotype_likelihood_cache.hpp"
#include "read_utils.hpp"

#include "haplotype_filter.hpp"

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
    
    using GenotypeLogLikelihood        = double;
    using SampleGenotypeLogLikelihoods = std::vector<GenotypeLogLikelihood>;
    using GenotypeLogLikelihoods       = std::unordered_map<SampleIdType, SampleGenotypeLogLikelihoods>;
    
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
    
    namespace debug {
        void print_haplotype_priors(const HaplotypePriorCounts& prior_counts, size_t n = 5);
        void print_haplotype_frequencies(const HaplotypeFrequencies& haplotype_frequencies, size_t n = 5);
        void print_genotype_log_likelihoods(const std::vector<Genotype<Haplotype>>& genotypes,
                                            const GenotypeLogLikelihoods& log_likelihoods, size_t n = 5);
        void print_genotype_log_marginals(const std::vector<Genotype<Haplotype>>& genotypes,
                                          const GenotypeLogMarginals& genotype_log_marginals,
                                          size_t n = 10);
        void print_genotype_posteriors(const std::vector<Genotype<Haplotype>>& genotypes,
                                       const GenotypePosteriors& genotype_posteriors,
                                       size_t n = 10);
        
        void print_read_haplotype_liklihoods(const std::vector<Haplotype>& haplotypes, const ReadMap& reads,
                                             HaplotypeLikelihoodCache& haplotype_likelihoods, size_t n = 3);
        void print_read_genotype_liklihoods(const std::vector<Genotype<Haplotype>>& genotypes, const ReadMap& reads,
                                 ReadModel& read_model, size_t n = 3);
    } // namespace debug
    
    // non member methods
    
    GenotypeLogLikelihoods
    compute_genotype_log_likelihoods(const std::vector<Genotype<Haplotype>>& genotypes,
                                     const ReadMap& reads, ReadModel& read_model)
    {
        using std::cbegin; using std::cend; using std::begin; using std::transform;
        
        GenotypeLogLikelihoods result {};
        result.reserve(reads.size());
        
        for (const auto& sample_reads_p : reads) {
            const auto& sample       = sample_reads_p.first;
            const auto& sample_reads = sample_reads_p.second;
            
            SampleGenotypeLogLikelihoods sample_log_likelihoods(genotypes.size());
            
            transform(cbegin(genotypes), cend(genotypes), begin(sample_log_likelihoods),
                      [&sample_reads, &read_model] (const auto& genotype) {
                          return read_model.log_probability(cbegin(sample_reads), cend(sample_reads), genotype);
                      });
            
            result.emplace(sample, std::move(sample_log_likelihoods));
        }
        
        return result;
    }
    
    GenotypeLogMarginals
    init_genotype_log_marginals(const std::vector<Genotype<Haplotype>>& genotypes,
                                const HaplotypeFrequencies& haplotype_frequencies)
    {
        GenotypeLogMarginals result {};
        result.reserve(genotypes.size());
        
        for (const auto& genotype : genotypes) {
            result.emplace_back(genotype, log_hardy_weinberg(genotype, haplotype_frequencies));
        }
        
        return result;
    }
    
    void update_genotype_log_marginals(GenotypeLogMarginals& current_log_marginals,
                                       const HaplotypeFrequencies& haplotype_frequencies)
    {
        using std::begin; using std::end; using std::for_each;
        
        for_each(begin(current_log_marginals), end(current_log_marginals),
                 [&haplotype_frequencies] (auto& p) {
                     p.log_probability = log_hardy_weinberg(p.genotype, haplotype_frequencies);
                 });
    }
    
    static void normalise_exp(SampleGenotypePosteriors& unnormalised_log_probabilities)
    {
        using std::cbegin; using std::cend; using std::begin; using std::transform;
        
        const auto norm = Maths::log_sum_exp<double>(unnormalised_log_probabilities);
        transform(cbegin(unnormalised_log_probabilities), cend(unnormalised_log_probabilities),
                  begin(unnormalised_log_probabilities), [norm] (double p) {
                      return std::exp(p - norm);
                  });
    }
    
    GenotypePosteriors
    init_genotype_posteriors(const GenotypeLogMarginals& genotype_log_marginals,
                             const GenotypeLogLikelihoods& genotype_log_likilhoods)
    {
        using std::cbegin; using std::cend; using std::begin; using std::transform;
        
        GenotypePosteriors result {};
        result.reserve(genotype_log_likilhoods.size());
        
        for (const auto& sample_genotype_log_likilhoods : genotype_log_likilhoods) {
            SampleGenotypePosteriors sample_result(genotype_log_marginals.size());
            
            transform(cbegin(genotype_log_marginals), cend(genotype_log_marginals),
                      cbegin(sample_genotype_log_likilhoods.second), begin(sample_result),
                      [] (const auto& genotype_log_marginal, const auto genotype_log_likilhood) {
                          return genotype_log_marginal.log_probability + genotype_log_likilhood;
                      });
            
            normalise_exp(sample_result);
            
            result.emplace(sample_genotype_log_likilhoods.first, std::move(sample_result));
        }
        
        return result;
    }
    
    void update_genotype_posteriors(GenotypePosteriors& current_genotype_posteriors,
                                    const HaplotypeFrequencies& haplotype_frequencies,
                                    const GenotypeLogMarginals& genotype_log_marginals,
                                    const GenotypeLogLikelihoods& genotype_log_likilhoods)
    {
        using std::cbegin; using std::cend; using std::begin; using std::transform;
        
        for (auto& sample_genotype_posteriors : current_genotype_posteriors) {
            
            const auto& sample_genotype_log_likilhoods = genotype_log_likilhoods.at(sample_genotype_posteriors.first);
            
            transform(cbegin(genotype_log_marginals), cend(genotype_log_marginals),
                      cbegin(sample_genotype_log_likilhoods), begin(sample_genotype_posteriors.second),
                      [] (const auto& genotype_log_marginal, const auto& genotype_log_likilhood) {
                          return genotype_log_marginal.log_probability + genotype_log_likilhood;
                      });
            
            normalise_exp(sample_genotype_posteriors.second);
        }
    }
    
    HaplotypeFrequencies
    init_haplotype_frequencies(const HaplotypePriorCounts& haplotype_prior_counts, const double prior_count_sum)
    {
        HaplotypeFrequencies result {};
        result.reserve(haplotype_prior_counts.size());
        
        for (const auto& haplotype_count : haplotype_prior_counts) {
            result.emplace(haplotype_count.first, haplotype_count.second / prior_count_sum);
        }
        
        return result;
    }
    
    double update_haplotype_frequencies(HaplotypeFrequencies& current_haplotype_frequencies,
                                        const HaplotypePriorCounts& haplotype_prior_counts,
                                        const std::vector<Genotype<Haplotype>>& genotypes,
                                        const GenotypePosteriors& genotype_posteriors,
                                        const double prior_count_sum)
    {
        using std::cbegin; using std::cend; using std::accumulate;
        
        double max_frequency_change {};
        
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
                           HaplotypeFrequencies& haplotype_frequencies,
                           GenotypeLogMarginals& genotype_log_marginals,
                           const GenotypeLogLikelihoods& genotype_log_likilhoods,
                           const HaplotypePriorCounts& haplotype_prior_counts,
                           const double prior_count_sum)
    {
        auto max_change = update_haplotype_frequencies(haplotype_frequencies, haplotype_prior_counts,
                                                       genotypes, genotype_posteriors, prior_count_sum);
        
        //debug::print_haplotype_frequencies(haplotype_frequencies);
        
        update_genotype_log_marginals(genotype_log_marginals, haplotype_frequencies);
        
        update_genotype_posteriors(genotype_posteriors, haplotype_frequencies,
                                   genotype_log_marginals, genotype_log_likilhoods);
        
        //debug::print_genotype_posteriors(genotypes, genotype_posteriors);
        
        return max_change;
    }
    
    Population::Latents
    make_latents(std::vector<Genotype<Haplotype>>&& genotypes, GenotypePosteriors&& genotype_posteriors,
                 HaplotypeFrequencies&& haplotype_frequencies)
    {
        Population::GenotypeProbabilities result {};
        result.reserve(genotype_posteriors.size()); // num samples
        
        for (auto&& sample_genotype_posteriors : genotype_posteriors) {
            Population::SampleGenotypeProbabilities sample_result {};
            sample_result.reserve(genotypes.size());
            
            auto it = std::cbegin(sample_genotype_posteriors.second);
            
            for (const auto& genotype : genotypes) {
                sample_result.emplace(genotype, *it);
                ++it;
            }
            
            result.emplace(std::move(sample_genotype_posteriors.first), std::move(sample_result));
        }
        
        return Population::Latents {std::move(result), std::move(haplotype_frequencies)};
    }
    
    // private methods
    
    Population::Latents
    Population::evaluate(const std::vector<Haplotype>& haplotypes, const ReadMap& reads, ReferenceGenome& reference)
    {
        if (haplotypes.size() == 1) {
            // TODO: catch this case to avoid computing
        }
        
//        auto filtered_haplotypes = filter_haplotypes(haplotypes, reads, 10);
//        
//        std::cout << "filtered_haplotypes  " << filtered_haplotypes.size() << std::endl;
//        for (const auto& h : filtered_haplotypes) {
//            print_variant_alleles(h);
//            std::cout << std::endl;
//        }
        
        HaplotypeLikelihoodCache haplotype_likelihoods {reads, haplotypes};
        
        auto genotypes = generate_all_genotypes(haplotypes, ploidy_);
        
        std::cout << "there are " << genotypes.size() << " candidate genotypes" << std::endl;
        
        ReadModel read_model {ploidy_, haplotype_likelihoods};
        
        //auto start = std::chrono::system_clock::now();
        
        const auto genotype_log_likilhoods = compute_genotype_log_likelihoods(genotypes, reads, read_model);
        
        //auto end = std::chrono::system_clock::now();
        //std::cout << std::chrono::duration_cast<std::chrono::milliseconds>(end - start).count() << "ms" << std::endl;
        
//        debug::print_genotype_log_likelihoods(genotypes, genotype_log_likilhoods, 20);
//        std::cout << std::endl;
//        debug::print_read_haplotype_liklihoods(haplotypes, reads, haplotype_likelihoods);
//        std::cout << std::endl;
//        debug::print_read_genotype_liklihoods(genotypes, reads, read_model, 10);
        
        haplotype_likelihoods.clear();
        
        auto haplotype_prior_counts = compute_haplotype_prior_counts(haplotypes, reference, haplotype_prior_model_);
        
        const auto prior_count_sum = Maths::sum_values(haplotype_prior_counts);
        
        auto haplotype_frequencies  = init_haplotype_frequencies(haplotype_prior_counts, prior_count_sum);
        auto genotype_log_marginals = init_genotype_log_marginals(genotypes, haplotype_frequencies);
        
//        debug::print_haplotype_priors(haplotype_prior_counts);
//        debug::print_haplotype_frequencies(haplotype_frequencies);
//        debug::print_genotype_log_marginals(genotypes, genotype_log_marginals);
        //exit(0);
        
        auto genotype_posteriors = init_genotype_posteriors(genotype_log_marginals, genotype_log_likilhoods);
        
        //debug::print_genotype_posteriors(genotypes, genotype_posteriors);
        //exit(0);
        
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
        
        void print_haplotype_priors(const HaplotypePriorCounts& prior_counts, size_t n)
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
        
        void print_haplotype_frequencies(const HaplotypeFrequencies& haplotype_frequencies, size_t n)
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
                                            const GenotypeLogLikelihoods& log_likelihoods, size_t n)
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
        
        void print_read_haplotype_liklihoods(const std::vector<Haplotype>& haplotypes, const ReadMap& reads,
                                             HaplotypeLikelihoodCache& haplotype_likelihoods, size_t n)
        {
            auto m = std::min(n, haplotypes.size());
            
            std::cout << "top " << m << " haplotype likelihoods for each read in each sample" << std::endl;
            
            for (const auto& sample_reads : reads) {
                std::cout << sample_reads.first << ":" << std::endl;
                
                for (const auto& read : sample_reads.second) {
                    std::cout << read.get_region() << " " << read.get_cigar_string() << ":" << std::endl;
                    
                    std::vector<std::pair<Haplotype, double>> top {};
                    top.reserve(haplotypes.size());
                    
                    for (const auto& haplotype : haplotypes) {
                        top.emplace_back(haplotype, haplotype_likelihoods.log_probability(read, haplotype));
                    }
                    
                    std::sort(std::begin(top), std::end(top), IsBigger<Haplotype, double>());
                    
                    for (unsigned i {}; i < m; ++i) {
                        std::cout << "\t* ";
                        print_variant_alleles(top[i].first);
                        std::cout << " " << std::setprecision(10) << top[i].second << std::endl;
                    }
                }
            }
        }
        
        void print_read_genotype_liklihoods(const std::vector<Genotype<Haplotype>>& genotypes, const ReadMap& reads,
                                 ReadModel& read_model, size_t n)
        {
            auto m = std::min(n, genotypes.size());
            
            std::cout << "top " << n << " genotype likelihoods for each read in each sample" << std::endl;
            
            for (const auto& sample_reads : reads) {
                std::cout << sample_reads.first << ":" << std::endl;
                for (const auto& read : sample_reads.second) {
                    std::cout << read.get_region() << " " << read.get_cigar_string() << ":" << std::endl;
                    
                    std::vector<std::pair<Genotype<Haplotype>, double>> top {};
                    top.reserve(genotypes.size());
                    
                    for (const auto& genotype : genotypes) {
                        top.emplace_back(genotype, read_model.log_probability(read, genotype));
                    }
                    std::sort(std::begin(top), std::end(top), IsBigger<Genotype<Haplotype>, double>());
                    
                    for (unsigned i {}; i < m; ++i) {
                        std::cout << "\t* ";
                        print_variant_alleles(top[i].first);
                        std::cout << " " << std::setprecision(10) << top[i].second << std::endl;
                    }
                }
            }
        }
        
    } // namespace debug
    } // namespace GenotypeModel
} // namespace Octopus

