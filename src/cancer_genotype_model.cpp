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
#include <cmath>

#include "single_read_model.hpp"
#include "common.hpp"
#include "maths.hpp"

#include <iostream> // TEST

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
    
    using HaplotypeFrequencies        = std::unordered_map<Haplotype, double>;
    using GenotypeMarginals           = std::unordered_map<Genotype<Haplotype>, double>;
    using SampleGenotypeWeightsCounts = std::array<double, 3>;
    using SampleGenotypeWeights       = std::array<double, 3>;
    using GenotypeWeightsCounts       = std::unordered_map<Octopus::SampleIdType, SampleGenotypeWeightsCounts>;
    using GenotypeWeights             = std::unordered_map<Octopus::SampleIdType, SampleGenotypeWeights>;
    
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
                                   const MappableSet<AlignedRead>& reads)
    {
        static SingleReadModel rm {1000};
        
        double result {};
        
        for (const auto read : reads) {
            result += Maths::log_sum_exp(std::log(genotype_weights[0]) + rm.log_probability(read, genotype[0]),
                                  std::log(genotype_weights[1]) + rm.log_probability(read, genotype[1]),
                                  std::log(genotype_weights[2]) + rm.log_probability(read, genotype[2]));
        }
        
        return result;
    }
    
    void normalise_exp(Cancer::GenotypeProbabilities& unnormalised_log_genotype_probabilities)
    {
        std::vector<double> tmp(unnormalised_log_genotype_probabilities.size());
        
        std::transform(std::cbegin(unnormalised_log_genotype_probabilities),
                       std::cend(unnormalised_log_genotype_probabilities), std::begin(tmp),
                       [] (const auto& p) { return p.second; });
        
        auto norm = Maths::log_sum_exp<double>(std::cbegin(tmp), std::cend(tmp));
        
        for (auto& p : unnormalised_log_genotype_probabilities) {
            p.second -= norm;
            p.second = std::exp(p.second);
        }
    }
    
    Cancer::GenotypeProbabilities
    compute_genotype_posteriors(const std::vector<CancerGenotype<Haplotype>>& genotypes, const ReadMap& reads,
                                const HaplotypeFrequencies& haplotype_frequencies,
                                const Cancer::GenotypeWeights& genotype_weights)
    {
        Cancer::GenotypeProbabilities result {};
        result.reserve(genotypes.size());
        
        for (const auto& genotype : genotypes) {
            auto p = log_hardy_weinberg(genotype, haplotype_frequencies);
            
            for (const auto& sample_reads : reads) {
                const auto& w = genotype_weights.at(sample_reads.first);
                p += genotype_log_likelihood(genotype, w, sample_reads.second);
                p -= sample_reads.second.size() * std::log(sum(w));
            }
            
            result.emplace(genotype, p);
        }
        
        normalise_exp(result);
        
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
    
    Cancer::GenotypeWeightResponsibilities
    compute_genotype_weight_responsibilities(const Cancer::GenotypeProbabilities& genotype_probabilities,
                                             const ReadMap& reads, const Cancer::GenotypeWeights& genotype_weights)
    {
        static SingleReadModel rm {1000};
        
        Cancer::GenotypeWeightResponsibilities result {};
        result.reserve(reads.size());
        
        for (const auto& sample_reads : reads) {
            std::vector<std::array<double, 3>> v {};
            v.reserve(sample_reads.second.size());
            
            auto log_weight = log(genotype_weights.at(sample_reads.first));
            
            for (const auto& read : sample_reads.second) {
                std::array<double, 3> p {0.0, 0.0, 0.0};
                
                for (unsigned k {}; k < 3; ++k) {
                    std::vector<double> lg {};
                    lg.reserve(genotype_probabilities.size());
                    
                    for (const auto& p : genotype_probabilities) {
                        lg.push_back(std::log(p.second) + log_weight[k] + rm.log_probability(read, p.first[k]));
                    }
                    
                    p[k] = Maths::log_sum_exp<double>(std::cbegin(lg), std::cend(lg));
                }
                
                normalise_exp(p);
                
                v.push_back(p);
            }
            
            result.emplace(sample_reads.first, std::move(v));
        }
        
        return result;
    }
    
    HaplotypeFrequencies
    compute_haplotype_frequencies(const std::vector<Haplotype>& haplotypes,
                                  const Cancer::GenotypeProbabilities& genotype_posteriors)
    {
        HaplotypeFrequencies result {};
        result.reserve(haplotypes.size());
        
        for (const auto& haplotype : haplotypes) {
            double p {};
            
            for (const auto& genotype_posterior : genotype_posteriors) {
                p += genotype_posterior.first.num_occurences(haplotype) * genotype_posterior.second;
            }
            
            result.emplace(haplotype, p / 3);
        }
        
        return result;
    }
    
    Cancer::GenotypeWeights
    compute_genotype_weights(const Cancer::GenotypeWeightsPriors& prior_counts,
                             const Cancer::GenotypeWeightResponsibilities& genotype_weight_responsibilities)
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
    
    double do_em_iteration(const std::vector<Haplotype>& haplotypes,
                           const std::vector<CancerGenotype<Haplotype>>& genotypes,
                           HaplotypeFrequencies& haplotype_frequencies,
                           const Cancer::GenotypeWeightsPriors& weight_priors,
                           Cancer::GenotypeWeights& genotype_weights,
                           Cancer::GenotypeProbabilities& genotypes_probabilities,
                           Cancer::GenotypeWeightResponsibilities& genotype_weight_responsibilities,
                           const ReadMap& reads)
    {
        auto new_haplotype_frequencies = compute_haplotype_frequencies(haplotypes, genotypes_probabilities);
        auto new_genotype_weights      = compute_genotype_weights(weight_priors, genotype_weight_responsibilities);
        
        genotypes_probabilities = compute_genotype_posteriors(genotypes, reads, new_haplotype_frequencies, genotype_weights);
        genotype_weight_responsibilities = compute_genotype_weight_responsibilities(genotypes_probabilities, reads, new_genotype_weights);
        
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
        auto haplotype_prior_counts = compute_haplotype_prior_counts(haplotypes, reference, haplotype_prior_model_);
        
        auto genotypes = generate_all_cancer_genotypes(haplotypes, 2);
        
        GenotypeWeightsPriors weight_priors {};
        for (const auto& s : reads) {
            if (s.first == normal_sample_) {
                weight_priors.emplace(s.first, SampleGenotypeWeightsPriors {1.0, 1.0, 0.01});
            } else {
                weight_priors.emplace(s.first, SampleGenotypeWeightsPriors {0.5, 0.5, 0.5});
            }
        }
        
        auto haplotype_frequencies = init_haplotype_frequencies(haplotype_prior_counts);
        auto genotype_weights      = init_genotype_weights(weight_priors);
        
        auto genotype_posteriors = compute_genotype_posteriors(genotypes, reads, haplotype_frequencies, genotype_weights);
        auto genotype_weight_responsibilities = compute_genotype_weight_responsibilities(genotype_posteriors, reads, genotype_weights);
        
        for (unsigned n {0}; n < max_em_iterations_; ++n) {
            //std::cout << "EM iteration " << n << std::endl;
            auto c = do_em_iteration(haplotypes, genotypes, haplotype_frequencies, weight_priors,
                                     genotype_weights, genotype_posteriors,
                                     genotype_weight_responsibilities, reads);
            if (c < em_epsilon_) break;
        }
        
        Latents result {};
        result.genotype_posteriors = std::move(genotype_posteriors);
        result.genotype_weights    = std::move(genotype_weights);
        
        return result;
    }
    
    } // namespace GenotypeModel
} // namespace Octopus
