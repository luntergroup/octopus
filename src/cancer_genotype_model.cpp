//
//  cancer_genotype_model.cpp
//  Octopus
//
//  Created by Daniel Cooke on 26/08/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "cancer_genotype_model.hpp"

#include <array>
#include <numeric>    // std::accumulate
#include <algorithm>  // std::transform, std::nth_element
#include <iterator>
#include <cmath>

#include "common.hpp"
#include "maths.hpp"
#include "read_utils.hpp"

#include "threaded_transform.hpp" // TEST

#include <iostream>   // DEBUG

namespace Octopus
{
    namespace GenotypeModel
    {
    // public methods
    
    Cancer::Cancer(SampleIdType normal_sample, unsigned max_em_iterations, double em_epsilon)
    :
    max_em_iterations_ {max_em_iterations},
    em_epsilon_ {em_epsilon},
    normal_sample_ {std::move(normal_sample)}
    {}
    
    // non member methods
    
    using HaplotypeFrequencies           = std::unordered_map<Haplotype, double>;
    using SampleGenotypeWeightsCounts    = std::array<double, 3>;
    using SampleGenotypeWeights          = std::array<double, 3>;
    using GenotypeWeightsCounts          = std::unordered_map<SampleIdType, SampleGenotypeWeightsCounts>;
    using GenotypeWeights                = std::unordered_map<SampleIdType, SampleGenotypeWeights>;
    using GenotypeWeightResponsibilities = std::unordered_map<SampleIdType, std::vector<std::array<double, 3>>>;
    using GenotypeLogPosterior           = double;
    using GenotypeLogPosteriors          = std::vector<GenotypeLogPosterior>;
    
    namespace debug {
        std::ostream& operator<<(std::ostream& os, const std::array<double, 3>& arr);
        std::ostream& operator<<(std::ostream& os, const std::unordered_map<Octopus::SampleIdType, std::array<double, 3>>& m);
        void print_top_genotypes(const std::vector<CancerGenotype<Haplotype>>& genotypes,
                                 const GenotypeLogPosteriors& genotype_log_posteriors,
                                 const size_t n = 20);
        void print_weight_responsabilities(const GenotypeWeightResponsibilities& responsabilities,
                                           const ReadMap& reads);
    } // namespace debug
    
    double sum(const std::array<double, 3>& arr)
    {
        return arr[0] + arr[1] + arr[2];
    }
    
    SampleGenotypeWeights expected_value(const SampleGenotypeWeightsCounts& counts)
    {
        auto n = sum(counts);
        return SampleGenotypeWeights {counts[0] / n, counts[1] / n, counts[2] / n};
    }
    
    Cancer::GenotypeWeights init_genotype_weights(const Cancer::GenotypeWeightsPriors& weight_counts)
    {
        Cancer::GenotypeWeights result {};
        result.reserve(weight_counts.size());
        
        for (const auto& sample_weight_counts : weight_counts) {
            const auto& curr = sample_weight_counts.second;
            const auto n = sum(curr);
            result.emplace(sample_weight_counts.first, Cancer::SampleGenotypeWeights {curr[0] / n, curr[1] / n, curr[2] / n});
        }
        
        return result;
    }
    
    std::array<double, 3> log(const std::array<double, 3>& arr)
    {
        return {std::log(arr[0]), std::log(arr[1]), std::log(arr[2])};
    }
    
    double genotype_log_probability(const CancerGenotype<Haplotype>& genotype,
                                    const HaplotypeFrequencies& haplotype_frequencies)
    {
        return log_hardy_weinberg(genotype.get_germline_genotype(), haplotype_frequencies) +
                std::log(haplotype_frequencies.at(genotype.get_cancer_element()));
    }
    
    double genotype_log_likelihood(const CancerGenotype<Haplotype>& genotype,
                                   const Cancer::SampleGenotypeWeights& genotype_weights,
                                   const MappableSet<AlignedRead>& reads, SingleReadModel& rm)
    {
        using std::cbegin; using std::cend; using std::accumulate;
        
        const auto log_weights = log(genotype_weights);
        
//        for (const auto& read : reads) {
//            std::cout << get_begin(read) << " " << read.get_cigar_string() << std::endl;
//            std::cout << "\t* " << rm.log_probability(read, genotype[0]) << std::endl;
//            std::cout << "\t* " << rm.log_probability(read, genotype[1]) << std::endl;
//            std::cout << "\t* " << rm.log_probability(read, genotype[2]) << std::endl;
//        }
        
        return accumulate(cbegin(reads), cend(reads), 0.0,
                          [&log_weights, &genotype, &rm] (double curr, const auto& read) {
                              return curr + Maths::log_sum_exp(log_weights[0] + rm.log_probability(read, genotype[0]),
                                                               log_weights[1] + rm.log_probability(read, genotype[1]),
                                                               log_weights[2] + rm.log_probability(read, genotype[2]));
                          });
    }
    
    GenotypeLogPosteriors
    init_genotype_log_posteriors(const std::vector<CancerGenotype<Haplotype>>& genotypes,
                                 const ReadMap& reads,
                                 const HaplotypeFrequencies& haplotype_frequencies,
                                 const Cancer::GenotypeWeights& genotype_weights,
                                 SingleReadModel& read_model)
    {
        using std::cbegin; using std::cend; using std::begin; using std::transform;
        
        GenotypeLogPosteriors result(genotypes.size());
        
        transform(cbegin(genotypes), cend(genotypes), begin(result),
                  [&haplotype_frequencies, &genotype_weights, &reads, &read_model]
                  (const auto& genotype) {
                      double p {genotype_log_probability(genotype, haplotype_frequencies)};
                      
                      for (const auto& sample_reads : reads) {
                          p += genotype_log_likelihood(genotype, genotype_weights.at(sample_reads.first),
                                                       sample_reads.second, read_model);
                      }
                      
                      return p;
                  });
        
        const auto norm = Maths::log_sum_exp<double>(result);
        
        for (auto& p : result) p -= norm;
        
        return result;
    }
    
    void remove_genotypes(std::vector<CancerGenotype<Haplotype>>& genotypes,
                          GenotypeLogPosteriors& genotype_log_posteriors,
                          size_t max_genotypes)
    {
        if (genotypes.empty() || genotypes.size() <= max_genotypes) return;
        
        auto tmp = genotype_log_posteriors;
        
        std::nth_element(std::begin(tmp), std::begin(tmp) + max_genotypes, std::end(tmp), std::greater<double>());
        
        double cutoff = tmp[max_genotypes];
        
        auto gitr = std::cbegin(genotypes);
        auto gend = std::cend(genotypes);
        auto pitr = std::cbegin(genotype_log_posteriors);
        
        std::vector<CancerGenotype<Haplotype>> kept_genotypes {};
        kept_genotypes.reserve(max_genotypes);
        GenotypeLogPosteriors kept_posteriors {};
        kept_posteriors.reserve(max_genotypes);
        
        for (; gitr != gend; ++gitr, ++pitr) {
            if (*pitr > cutoff) {
                kept_genotypes.push_back(*gitr);
                kept_posteriors.push_back(*pitr);
            }
        }
        
        genotypes = kept_genotypes;
        genotype_log_posteriors = kept_posteriors;
    }
    
    void update_genotype_log_posteriors(GenotypeLogPosteriors& current,
                                        const std::vector<CancerGenotype<Haplotype>>& genotypes,
                                        const ReadMap& reads,
                                        const HaplotypeFrequencies& haplotype_frequencies,
                                        const Cancer::GenotypeWeights& genotype_weights,
                                        SingleReadModel& read_model)
    {
        using std::cbegin; using std::cend; using std::begin; using std::transform;
        
        std::vector<double> ps {};
        ps.reserve(current.size());
        
        transform(cbegin(genotypes), cend(genotypes), cbegin(current), begin(current),
                  [&haplotype_frequencies, &genotype_weights, &reads, &read_model]
                  (const auto& genotype, double curr) {
                      double p {genotype_log_probability(genotype, haplotype_frequencies)};
                      
                      for (const auto& sample_reads : reads) {
                          p += genotype_log_likelihood(genotype, genotype_weights.at(sample_reads.first),
                                                       sample_reads.second, read_model);
                      }
                      
                      return p;
                  });
        
        const auto norm = Maths::log_sum_exp<double>(current);
        
        for (auto& p : current) p -= norm;
    }
    
    void normalise_exp(std::array<double, 3>& arr)
    {
        auto norm = Maths::log_sum_exp(arr[0], arr[1], arr[2]);
        
        for (auto& e : arr) {
            e -= norm;
            e = std::exp(e);
        }
    }
    
    GenotypeWeightResponsibilities
    init_genotype_weight_responsibilities(const GenotypeLogPosteriors& genotype_posteriors,
                                          const std::vector<CancerGenotype<Haplotype>>& genotypes,
                                          const Cancer::GenotypeWeights& genotype_weights,
                                          const ReadMap& reads, SingleReadModel& rm)
    {
        using std::cbegin; using std::cend; using std::begin; using std::transform;
        
        GenotypeWeightResponsibilities result {};
        result.reserve(reads.size());
        
        for (const auto& sample_reads : reads) {
            std::vector<std::array<double, 3>> v {};
            v.reserve(sample_reads.second.size());
            
            auto log_weights = log(genotype_weights.at(sample_reads.first));
            
            for (const auto& read : sample_reads.second) {
                std::array<double, 3> p {0.0, 0.0, 0.0};
                
                for (unsigned k {}; k < 3; ++k) {
                    std::vector<double> lg(genotypes.size());
                    
                    transform(cbegin(genotypes), cend(genotypes), cbegin(genotype_posteriors),
                              begin(lg), [&read, &log_weights, &rm, k]
                              (const auto& genotype, double log_posterior) {
                                  return log_posterior + log_weights[k] + rm.log_probability(read, genotype[k]);
                              });
                    
                    p[k] = Maths::log_sum_exp<double>(lg);
                }
                
                normalise_exp(p);
                
                v.push_back(p);
            }
            
            result.emplace(sample_reads.first, std::move(v));
        }
        
        return result;
    }
    
    void update_genotype_weight_responsibilities(GenotypeWeightResponsibilities& current,
                                                 const std::vector<CancerGenotype<Haplotype>>& genotypes,
                                                 const GenotypeLogPosteriors& genotype_log_posteriors,
                                                 const Cancer::GenotypeWeights& genotype_weights,
                                                 const ReadMap& reads,
                                                 SingleReadModel& rm)
    {
        using std::cbegin; using std::cend; using std::begin; using std::transform;
        
        for (auto& sample_p : current) {
            const auto& sample = sample_p.first;
            
            auto log_weights = log(genotype_weights.at(sample));
            
            auto itr = std::begin(sample_p.second);
            
            for (const auto& read : reads.at(sample)) {
                std::array<double, 3> p {0.0, 0.0, 0.0};
                
                for (unsigned k {}; k < 3; ++k) {
                    std::vector<double> lg(genotypes.size());
                    
                    transform(cbegin(genotypes), cend(genotypes), cbegin(genotype_log_posteriors),
                              begin(lg), [&read, &log_weights, &rm, k]
                              (const auto& genotype, double log_posterior) {
                                  return log_posterior + log_weights[k] + rm.log_probability(read, genotype[k]);
                              });
                    
                    p[k] = Maths::log_sum_exp<double>(lg);
                }
                
                normalise_exp(p);
                
                *itr = p;
                ++itr;
            }
        }
    }
    
    double update_haplotype_frequencies(HaplotypeFrequencies& current,
                                        const HaplotypePriorCounts& haplotype_prior_counts,
                                        const GenotypeLogPosteriors& genotype_log_posteriors,
                                        const std::vector<CancerGenotype<Haplotype>>& genotypes)
    {
        double max_frequency_change {0.0};
        
        const auto norm = 3 + Maths::sum_values(haplotype_prior_counts);
        
        for (auto& haplotype_frequency : current) {
            double p {haplotype_prior_counts.at(haplotype_frequency.first)};
            
            auto itr = std::cbegin(genotype_log_posteriors);
            
            for (const auto& genotype : genotypes) {
                p += genotype.count(haplotype_frequency.first) * std::exp(*itr);
                ++itr;
            }
            
            const auto new_frequency = p / norm;
            
            const auto curr_fequency_change = std::abs(haplotype_frequency.second - new_frequency);
            
            if (curr_fequency_change > max_frequency_change) max_frequency_change = curr_fequency_change;
            
            haplotype_frequency.second = new_frequency;
        }
        
        return max_frequency_change;
    }
    
    Cancer::GenotypeWeights
    compute_genotype_weights(const Cancer::GenotypeWeightsPriors& prior_counts,
                             const GenotypeWeightResponsibilities& genotype_weight_responsibilities)
    {
        Cancer::GenotypeWeights result {};
        result.reserve(prior_counts.size());
        
        for (const auto& sample_prior_counts : prior_counts) {
            const auto& sample = sample_prior_counts.first;
            const auto& z = genotype_weight_responsibilities.at(sample);
            
            Cancer::SampleGenotypeWeights curr {};
            
            auto norm = z.size() + sum(sample_prior_counts.second);
            
            for (unsigned k {}; k < 3; ++k) {
                curr[k] = std::accumulate(std::cbegin(z), std::cend(z), 0.0,
                                          [k] (double v, const auto& a) { return v + a[k]; });
                curr[k] += sample_prior_counts.second[k];
                curr[k] /= norm;
            }
            
            result.emplace(sample, curr);
        }
        
        return result;
    }
    
    double update_genotype_weights(Cancer::GenotypeWeights& current,
                                 const Cancer::GenotypeWeightsPriors& prior_counts,
                                 const GenotypeWeightResponsibilities& genotype_weight_responsibilities)
    {
        double max_genotype_weight_change {0.0};
        
        for (auto& sample_weights : current) {
            const auto& sample = sample_weights.first;
            const auto& z = genotype_weight_responsibilities.at(sample);
            
            Cancer::SampleGenotypeWeights curr {};
            
            const auto& sample_prior_counts = prior_counts.at(sample);
            
            auto norm = z.size() + sum(sample_prior_counts);
            
            for (unsigned k {}; k < 3; ++k) {
                curr[k] = std::accumulate(std::cbegin(z), std::cend(z), 0.0,
                                          [k] (double v, const auto& a) { return v + a[k]; });
                curr[k] += sample_prior_counts[k];
                curr[k] /= norm;
                
                const auto& curr_weight_change = std::abs(curr[k] - sample_weights.second[k]);
                
                if (curr_weight_change > max_genotype_weight_change) {
                    max_genotype_weight_change = curr_weight_change;
                }
            }
            
            sample_weights.second = curr;
        }
        
        return max_genotype_weight_change;
    }
    
    double do_em_iteration(const std::vector<CancerGenotype<Haplotype>>& genotypes,
                           HaplotypeFrequencies& haplotype_frequencies,
                           const Cancer::GenotypeWeightsPriors& weight_priors,
                           Cancer::GenotypeWeights& genotype_weights,
                           GenotypeLogPosteriors& genotypes_log_probabilities,
                           GenotypeWeightResponsibilities& genotype_weight_responsibilities,
                           const ReadMap& reads,
                           const HaplotypePriorCounts& haplotype_prior_counts,
                           SingleReadModel& read_model)
    {
        auto max_frequency_change = update_haplotype_frequencies(haplotype_frequencies, haplotype_prior_counts,
                                                                 genotypes_log_probabilities, genotypes);
        
        auto max_weight_change = update_genotype_weights(genotype_weights, weight_priors,
                                                         genotype_weight_responsibilities);
        
        update_genotype_log_posteriors(genotypes_log_probabilities, genotypes, reads,
                                       haplotype_frequencies, genotype_weights, read_model);
        
        update_genotype_weight_responsibilities(genotype_weight_responsibilities, genotypes,
                                                genotypes_log_probabilities, genotype_weights,
                                                reads, read_model);
        
        return std::max(max_frequency_change, max_weight_change);
    }
    
    // private methods
    
    Cancer::Latents
    Cancer::evaluate(const std::vector<Haplotype>& haplotypes, const ReadMap& reads, ReferenceGenome& reference)
    {
        read_model_ = SingleReadModel {max_sample_read_count(reads), haplotypes.size()};
        
        auto haplotype_prior_counts = compute_haplotype_prior_counts(haplotypes, reference, haplotype_prior_model_);
        
        auto genotypes = generate_all_cancer_genotypes(haplotypes, 2); // diploid only for now
        
        std::cout << "there are " << genotypes.size() << " candidate cancer genotypes" << std::endl;
        
        GenotypeWeightsPriors weight_priors {};
        for (const auto& s : reads) {
            if (s.first == normal_sample_) {
                weight_priors.emplace(s.first, SampleGenotypeWeightsPriors {1000.0, 1000.0, 0.0001});
            } else {
                weight_priors.emplace(s.first, SampleGenotypeWeightsPriors {1.0, 1.0, 0.8});
            }
        }
        
        auto haplotype_frequencies = init_haplotype_frequencies(haplotype_prior_counts);
        
        for (const auto& hf : haplotype_frequencies) {
            print_variant_alleles(hf.first);
            std::cout << " : " << hf.second << std::endl;
        }
        
        auto genotype_weights      = init_genotype_weights(weight_priors);
        
        std::cout << "prior mixture weights" << std::endl;
        for (const auto& w : genotype_weights) {
            std::cout << w.first << ": " << w.second[0] << " " << w.second[1] << " " << w.second[2] << std::endl;
        }
        
        auto genotype_log_posteriors = init_genotype_log_posteriors(genotypes, reads, haplotype_frequencies,
                                                                    genotype_weights, read_model_);
        
        remove_genotypes(genotypes, genotype_log_posteriors, 20);
        
        debug::print_top_genotypes(genotypes, genotype_log_posteriors);
        //exit(0);
        
        auto genotype_weight_responsibilities = init_genotype_weight_responsibilities(genotype_log_posteriors,
                                                                                      genotypes,
                                                                                      genotype_weights, reads,
                                                                                      read_model_);
        
        //debug::print_weight_responsabilities(genotype_weight_responsibilities, reads);
        
        for (unsigned n {0}; n < max_em_iterations_; ++n) {
            std::cout << "EM iteration " << n << std::endl;
            if (do_em_iteration(genotypes, haplotype_frequencies, weight_priors,
                                       genotype_weights, genotype_log_posteriors,
                                       genotype_weight_responsibilities, reads,
                                       haplotype_prior_counts, read_model_) < em_epsilon_) break;
        }
        
        Cancer::GenotypeProbabilities genotype_posteriors {};
        genotype_posteriors.reserve(genotypes.size());
        
        auto itr = std::cbegin(genotype_log_posteriors);
        
        for (auto&& genotype : genotypes) {
            genotype_posteriors.emplace(std::move(genotype), std::exp(*itr));
            ++itr;
        }
        
        return Latents {std::move(genotype_posteriors), std::move(genotype_weights)};
    }
    
    namespace debug {
        
        std::ostream& operator<<(std::ostream& os, const std::array<double, 3>& arr)
        {
            os << arr[0] << " " << arr[1] << " " << arr[2];
            return os;
        }
        
        std::ostream& operator<<(std::ostream& os, const std::unordered_map<Octopus::SampleIdType, std::array<double, 3>>& m)
        {
            for (const auto& p : m) os << p.first << ": " << p.second << "\n";
            return os;
        }
        
        void print_top_genotypes(const std::vector<CancerGenotype<Haplotype>>& genotypes,
                                 const GenotypeLogPosteriors& genotype_log_posteriors,
                                 const size_t n)
        
        {
            std::vector<std::pair<const CancerGenotype<Haplotype>*, double>> v {};
            v.reserve(genotypes.size());
            
            std::transform(std::cbegin(genotypes), std::cend(genotypes), std::cbegin(genotype_log_posteriors),
                           std::back_inserter(v),
                           [] (const auto& genotype, double glp) {
                               return std::make_pair(&genotype, glp);
                           });
            
            std::sort(std::begin(v), std::end(v), [] (const auto& lhs, const auto& rhs) {
                return lhs.second > rhs.second;
            });
            
            auto m = std::min(genotype_log_posteriors.size(), n);
            
            std::cout << "DEBUG: print top " << m << " log genotype posteriors" << std::endl;
            
            for (unsigned i {}; i < m; ++i) {
                print_variant_alleles(*v[i].first);
                std::cout << " " << v[i].second << std::endl;
            }
        }
        
        void print_weight_responsabilities(const GenotypeWeightResponsibilities& responsabilities,
                                           const ReadMap& reads)
        {
            std::cout << "DEBUG: printing all read responsabilities" << std::endl;
            
            for (const auto& sample_reads : reads) {
                std::cout << sample_reads.first << ": " << std::endl;
                
                auto read_itr     = std::cbegin(sample_reads.second);
                auto read_end_itr = std::cend(sample_reads.second);
                auto r_itr = std::cbegin(responsabilities.at(sample_reads.first));
                
                for (; read_itr != read_end_itr; ++read_itr, ++r_itr) {
                    std::cout << read_itr->get_region() << " " << read_itr->get_cigar_string() << " " << *r_itr << std::endl;
                }
            }
        }
    } // namespace debug
    
    } // namespace GenotypeModel
} // namespace Octopus
