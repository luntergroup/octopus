//
//  cancer_genotype_model.cpp
//  Octopus
//
//  Created by Daniel Cooke on 26/08/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "cancer_genotype_model.hpp"

#include <array>
#include <numeric>
#include <algorithm>
#include <iterator>
#include <cmath>

#include "common.hpp"
#include "maths.hpp"
#include "read_utils.hpp"

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
    
    struct GenotypeLogPosterior
    {
        GenotypeLogPosterior() = delete;
        GenotypeLogPosterior(const CancerGenotype<Haplotype>& genotype, double log_posterior)
        : genotype {genotype}, log_posterior {log_posterior} {}
        
        const CancerGenotype<Haplotype>& genotype;
        double log_posterior;
    };
    
    using GenotypeLogPosteriors = std::vector<GenotypeLogPosterior>;
    
        namespace Debug {
    
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
                                     const size_t n = 20)
            
            {
                std::vector<std::pair<const CancerGenotype<Haplotype>*, double>> v {};
                v.reserve(genotypes.size());
                
                std::transform(std::cbegin(genotypes), std::cend(genotypes),
                               std::cbegin(genotype_log_posteriors),
                               std::back_inserter(v),
                               [] (const auto& genotype, const auto& glp) {
                                   return std::make_pair(&genotype, glp.log_posterior);
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
        
        } // namespace DEBUG
    
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
    
    double genotype_log_likelihood(const CancerGenotype<Haplotype>& genotype,
                                   const Cancer::SampleGenotypeWeights& genotype_weights,
                                   const MappableSet<AlignedRead>& reads, SingleReadModel& rm)
    {
        double result {};
        
        for (const auto& read : reads) {
            result += Maths::log_sum_exp(std::log(genotype_weights[0]) + rm.log_probability(read, genotype[0]),
                                         std::log(genotype_weights[1]) + rm.log_probability(read, genotype[1]),
                                         std::log(genotype_weights[2]) + rm.log_probability(read, genotype[2]));
        }
        
        return result;
    }
    
    GenotypeLogPosteriors
    compute_genotype_log_posteriors(const std::vector<CancerGenotype<Haplotype>>& genotypes,
                                    const ReadMap& reads,
                                    const HaplotypeFrequencies& haplotype_frequencies,
                                    const Cancer::GenotypeWeights& genotype_weights,
                                    SingleReadModel& read_model)
    {
        GenotypeLogPosteriors result {};
        result.reserve(genotypes.size());
        
        std::vector<double> ps {};
        ps.reserve(genotypes.size());
        
        for (const auto& genotype : genotypes) {
            double p {log_hardy_weinberg(genotype, haplotype_frequencies)};
            
            for (const auto& sample_reads : reads) {
                p += genotype_log_likelihood(genotype, genotype_weights.at(sample_reads.first),
                                             sample_reads.second, read_model);
                p -= std::log(sample_reads.second.size());
            }
            
            result.emplace_back(genotype, p);
            ps.push_back(p);
        }
        
        auto norm = Maths::log_sum_exp<double>(ps);
        
        for (auto& p : result) p.log_posterior -= norm;
        
        return result;
    }
    
    std::array<double, 3> log(std::array<double, 3> arr)
    {
        for (auto& e : arr) e = std::log(e);
        return arr;
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
    compute_genotype_weight_responsibilities(const GenotypeLogPosteriors& genotype_log_probabilities,
                                             const Cancer::GenotypeWeights& genotype_weights,
                                             const ReadMap& reads,
                                             SingleReadModel& rm)
    {
        GenotypeWeightResponsibilities result {};
        result.reserve(reads.size());
        
        for (const auto& sample_reads : reads) {
            std::vector<std::array<double, 3>> v {};
            v.reserve(sample_reads.second.size());
            
            auto log_weights = log(genotype_weights.at(sample_reads.first));
            
            for (const auto& read : sample_reads.second) {
                std::array<double, 3> p {0.0, 0.0, 0.0};
                
                for (unsigned k {}; k < 3; ++k) {
                    std::vector<double> lg {};
                    lg.reserve(genotype_log_probabilities.size());
                    
                    for (const auto& p : genotype_log_probabilities) {
                        lg.push_back(p.log_posterior + log_weights[k] + rm.log_probability(read, p.genotype[k]));
                    }
                    
                    p[k] = Maths::log_sum_exp<double>(lg);
                }
                
                normalise_exp(p);
                
                v.push_back(p);
            }
            
            result.emplace(sample_reads.first, std::move(v));
        }
        
        return result;
    }
    
    HaplotypeFrequencies
    compute_haplotype_frequencies(const HaplotypePriorCounts& haplotype_prior_counts,
                                  const GenotypeLogPosteriors& genotype_log_probabilities)
    {
        HaplotypeFrequencies result {};
        result.reserve(haplotype_prior_counts.size());
        
        const auto norm = 3 + Maths::sum_values(haplotype_prior_counts);
        
        for (const auto& haplotype_count : haplotype_prior_counts) {
            double p {haplotype_count.second};
            
            for (const auto& glp : genotype_log_probabilities) {
                p += glp.genotype.count(haplotype_count.first) * std::exp(glp.log_posterior);
            }
            
            result.emplace(haplotype_count.first, p / norm);
        }
        
        return result;
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
    
    double max_genotype_weight_change(const Cancer::GenotypeWeights& old_weights,
                                      const Cancer::GenotypeWeights& new_weights)
    {
        double result {};
        
        for (const auto& w : new_weights) {
            for (unsigned k {}; k < 3; ++k) {
                auto change = std::abs(w.second[k] - old_weights.at(w.first)[k]);
                if (change > result) result = change;
            }
        }
        
        return result;
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
        auto new_haplotype_frequencies = compute_haplotype_frequencies(haplotype_prior_counts, genotypes_log_probabilities);
        auto new_genotype_weights      = compute_genotype_weights(weight_priors, genotype_weight_responsibilities);
        
        genotypes_log_probabilities = compute_genotype_log_posteriors(genotypes, reads, new_haplotype_frequencies,
                                                                      genotype_weights, read_model);
        genotype_weight_responsibilities = compute_genotype_weight_responsibilities(genotypes_log_probabilities, new_genotype_weights,
                                                                                    reads, read_model);
        
        auto c1 = max_haplotype_frequency_change(haplotype_frequencies, new_haplotype_frequencies);
        auto c2 = max_genotype_weight_change(genotype_weights, new_genotype_weights);
        
        haplotype_frequencies = new_haplotype_frequencies;
        genotype_weights      = new_genotype_weights;
        
        return std::max(c1, c2);
    }
    
    // private methods
    
    Cancer::Latents
    Cancer::evaluate(const std::vector<Haplotype>& haplotypes, const ReadMap& reads, ReferenceGenome& reference)
    {
        read_model_ = SingleReadModel {max_sample_read_count(reads), haplotypes.size()};
        
        auto haplotype_prior_counts = compute_haplotype_prior_counts(haplotypes, reference, haplotype_prior_model_);
        
        auto genotypes = generate_all_cancer_genotypes(haplotypes, 2);
        
        std::cout << "there are " << genotypes.size() << " candidate cancer genotypes" << std::endl;
        
        GenotypeWeightsPriors weight_priors {};
        for (const auto& s : reads) {
            if (s.first == normal_sample_) {
                weight_priors.emplace(s.first, SampleGenotypeWeightsPriors {100.0, 100.0, 0.01});
            } else {
                weight_priors.emplace(s.first, SampleGenotypeWeightsPriors {0.6, 0.6, 0.5});
            }
        }
        
        auto haplotype_frequencies = init_haplotype_frequencies(haplotype_prior_counts);
        auto genotype_weights      = init_genotype_weights(weight_priors);
        
        auto genotype_log_posteriors  = compute_genotype_log_posteriors(genotypes, reads, haplotype_frequencies,
                                                                        genotype_weights, read_model_);
        
        Debug::print_top_genotypes(genotypes, genotype_log_posteriors, 10);
        
        //exit(0);
        
        auto genotype_weight_responsibilities = compute_genotype_weight_responsibilities(genotype_log_posteriors,
                                                                                         genotype_weights, reads, read_model_);
        
        for (unsigned n {0}; n < max_em_iterations_; ++n) {
            std::cout << "EM iteration " << n << std::endl;
            double c = do_em_iteration(genotypes, haplotype_frequencies, weight_priors,
                                       genotype_weights, genotype_log_posteriors,
                                       genotype_weight_responsibilities, reads,
                                       haplotype_prior_counts, read_model_);
            if (c < em_epsilon_) break;
        }
        
        Debug::print_weight_responsabilities(genotype_weight_responsibilities, reads);
        
        Cancer::GenotypeProbabilities genotype_posteriors {};
        genotype_posteriors.reserve(genotypes.size());
        
        for (const auto& g : genotype_log_posteriors) {
            genotype_posteriors.emplace(g.genotype, std::exp(g.log_posterior));
        }
        
        return Latents {std::move(genotype_posteriors), std::move(genotype_weights)};
    }
    
    } // namespace GenotypeModel
} // namespace Octopus
