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
#include "read_model.hpp"
#include "common.hpp"
#include "maths.hpp"

#include <iostream> // TEST

namespace Octopus
{
    // public methods
    
    CancerGenotypeModel::CancerGenotypeModel(unsigned num_samples, SampleIdType normal_sample_id,
                                             unsigned max_em_iterations, double em_epsilon)
    :
    max_em_iterations_ {max_em_iterations},
    em_epsilon_ {em_epsilon},
    num_samples_ {num_samples},
    normal_sample_id_ {std::move(normal_sample_id)}
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
    
    static double log_hardy_weinberg(const Genotype<Haplotype>& genotype,
                                     const HaplotypeFrequencies& haplotype_frequencies)
    {
        auto unique_haplotypes = genotype.get_unique();
        
        std::vector<unsigned> occurences {};
        occurences.reserve(unique_haplotypes.size());
        
        double r {};
        
        for (const auto& haplotype : unique_haplotypes) {
            auto num_occurences = genotype.num_occurences(haplotype);
            occurences.push_back(num_occurences);
            r += num_occurences * std::log(haplotype_frequencies.at(haplotype));
        }
        
        return log_multinomial_coefficient<double>(occurences.cbegin(), occurences.cend()) + r;
    }
    
    static HaplotypeFrequencies init_haplotype_frequencies(const std::vector<Haplotype>& haplotypes)
    {
        HaplotypeFrequencies result {};
        result.reserve(haplotypes.size());
        
        const double uniform {1.0 / haplotypes.size()};
        
        for (const auto& haplotype : haplotypes) {
            result.emplace(haplotype, uniform);
        }
        
        return result;
    }
    
    GenotypeWeightsCounts
    init_genotype_weight_counts(const ReadMap& reads)
    {
        GenotypeWeightsCounts result {};
        
        for (const auto& sample_reads : reads) {
            result.emplace(sample_reads.first, SampleGenotypeWeightsCounts {0.5, 0.5, 0.5});
        }
        
        return result;
    }
    
    double genotype_log_likelihood(const Genotype<Haplotype>& genotype, const SampleGenotypeWeights& f,
                                   const MappableSet<AlignedRead>& reads)
    {
        static SingleReadModel rm {1000};
        
        double result {};
        
        for (const auto read : reads) {
            result += log_sum_exp(std::log(f[0]) + rm.log_probability(read, genotype[0]),
                                  std::log(f[1]) + rm.log_probability(read, genotype[1]),
                                  std::log(f[2]) + rm.log_probability(read, genotype[2]));
        }
        
        return result;
    }
    
    double log_joint_liklihood(const Genotype<Haplotype>& g, const GenotypeWeights& fs,
                               const ReadMap& reads, const HaplotypeFrequencies& pi,
                               const GenotypeWeightsCounts& alphas)
    {
        double p {};
        
        for (const auto sample_reads : reads) {
            auto f = fs.at(sample_reads.first);
            auto a = alphas.at(sample_reads.first);
            p += log_dirichlet<double>(a.begin(), a.end(), f.begin()) + genotype_log_likelihood(g, f, sample_reads.second);
        }
        
        return log_hardy_weinberg(g, pi) + p;
    }
    
    GenotypeMarginals
    compute_genotype_posteriors(const std::vector<Genotype<Haplotype>>& genotypes, const ReadMap& reads,
                                const HaplotypeFrequencies& pi, const GenotypeWeightsCounts& alphas)
    {
        GenotypeMarginals result {};
        result.reserve(genotypes.size());
        
        std::vector<double> t {};
        t.reserve(genotypes.size());
        
        for (const auto& g : genotypes) {
            auto p = log_hardy_weinberg(g, pi);
            
            for (const auto& sample_reads : reads) {
                const auto& s = sample_reads.first;
                const auto a  = alphas.at(s);
                p += genotype_log_likelihood(g, a, sample_reads.second);
                p -= sample_reads.second.size() * std::log(a[0] +  a[1] +  a[2]);
            }
            
            result.emplace(g, p);
            t.emplace_back(p);
        }
        
        auto norm = log_sum_exp<double>(t.cbegin(), t.cend());
        
        for (auto& g : result) {
            g.second -= norm;
            g.second = std::exp(g.second);
        }
        
        return result;
    }
    
    HaplotypeFrequencies
    compute_haplotype_frequencies(const GenotypeMarginals& genotype_posteriors)
    {
        HaplotypeFrequencies result {};
        
        for (const auto& posterior : genotype_posteriors) {
            for (auto haplotype : posterior.first.get_unique()) {
                result[haplotype] += posterior.first.num_occurences(haplotype) * posterior.second;
            }
        }
        
        for (auto& h : result) {
            h.second /= 3;
        }
        
        return result;
    }
    
    GenotypeWeightsCounts
    compute_component_pseudo_counts(const GenotypeWeightsCounts& prior_counts,
                                    const GenotypeMarginals& genotype_posteriors,
                                    const ReadMap& reads)
    {
        GenotypeWeightsCounts result {};
        result.reserve(prior_counts.size());
        
        SingleReadModel rm {1000};
        
        for (const auto& sample_prior_counts : prior_counts) {
            result[sample_prior_counts.first] = sample_prior_counts.second;
            
            for (const auto& genotype_posterior : genotype_posteriors) {
                double a {}, b {}, c {};
                
                for (const auto& read : reads.at(sample_prior_counts.first)) {
                    a += rm.log_probability(read, genotype_posterior.first.at(0));
                    b += rm.log_probability(read, genotype_posterior.first.at(1));
                    c += rm.log_probability(read, genotype_posterior.first.at(2));
                }
                
                auto n = log_sum_exp(a, b, c);
                a -= n;
                b -= n;
                c -= n;
                
                result[sample_prior_counts.first][0] += genotype_posterior.second * std::exp(a);
                result[sample_prior_counts.first][1] += genotype_posterior.second * std::exp(b);
                result[sample_prior_counts.first][2] += genotype_posterior.second * std::exp(c);
            }
        }
        
        return result;
    }
    
    void remove_redundant_genotypes(std::vector<Genotype<Haplotype>>& genotypes)
    {
        genotypes.erase(std::remove_if(std::begin(genotypes), std::end(genotypes),
                                       [] (const auto& genotype) {
                                           return genotype.is_homozygous();
                                       }), std::end(genotypes));
    }
    
    std::vector<std::array<unsigned, 3>> generate_all_ratios(unsigned n)
    {
        std::vector<std::array<unsigned, 3>> result {};
        
        for (unsigned i {}; i <= n; ++i) {
            for (unsigned j {i}; j <= n; ++j) {
                std::array<unsigned, 3> arr {i, j - i, n - j};
                if (std::find(std::cbegin(arr), std::cend(arr), n) == std::cend(arr)) {
                    result.push_back(arr);
                }
            }
        }
        
        return result;
    }
    
    std::vector<SampleGenotypeWeights> generate_all_fs(unsigned n)
    {
        auto ratios = generate_all_ratios(n);
        std::vector<SampleGenotypeWeights> result(ratios.size());
        double epsilon {std::nextafter(0, 1.0)};
        std::transform(std::cbegin(ratios), std::cend(ratios), std::begin(result),
                       [n, epsilon] (const auto& arr) {
            return SampleGenotypeWeights {
                static_cast<double>(arr[0]) / n + epsilon,
                static_cast<double>(arr[1]) / n + epsilon,
                static_cast<double>(arr[2]) / n + epsilon
            };
        });
        return result;
    }
    
    static double max_haplotype_frequency_change(const HaplotypeFrequencies& old_frequencies,
                                                 const HaplotypeFrequencies& new_frequencies)
    {
        double result {};
        
        for (const auto& h : new_frequencies) {
            auto change = std::abs(h.second - old_frequencies.at(h.first));
            if (change > result) result = change;
        }
        
        return result;
    }
    
    double do_em_iteration(GenotypeMarginals& genotype_posteriors,
                           HaplotypeFrequencies& haplotype_frequencies,
                           GenotypeWeightsCounts& weight_counts,
                           const ReadMap& reads)
    {
        auto new_haplotype_frequencies  = compute_haplotype_frequencies(genotype_posteriors);
        auto new_genotype_weight_counts = compute_component_pseudo_counts(weight_counts, genotype_posteriors, reads);
        
        auto max_change = max_haplotype_frequency_change(haplotype_frequencies, new_haplotype_frequencies);
        
        
        
        haplotype_frequencies = new_haplotype_frequencies;
        
        return max_change;
    }
    
    // private methods
    
    CancerGenotypeModel::GenotypeProbabilities
    CancerGenotypeModel::do_evaluate(const std::vector<Haplotype>& haplotypes, const ReadMap& reads)
    {
        auto ratios = generate_all_fs(20);
        
        GenotypeProbabilities result {};
        
        auto haplotype_frequencies  = init_haplotype_frequencies(haplotypes);
        auto genotype_weight_counts = init_genotype_weight_counts(reads);
        
        auto genotypes = generate_all_genotypes(haplotypes, 3);
        
        remove_redundant_genotypes(genotypes);
        
        auto genotype_posteriors = compute_genotype_posteriors(genotypes, reads, haplotype_frequencies, genotype_weight_counts);
        
        for (unsigned n {0}; n < max_em_iterations_; ++n) {
            std::cout << "EM iteration " << n << std::endl;
            auto c = do_em_iteration(genotype_posteriors, haplotype_frequencies, genotype_weight_counts, reads);
            if (c < em_epsilon_) break;
        }
        
        for (const auto& a : genotype_weight_counts) {
            std::cout << a.first << " " << expected_value(a.second) << std::endl;
        }
        
        for (const auto& h : haplotype_frequencies) {
            std::cout << h.first << " " << h.second << std::endl;
        }
        
        auto it = std::max_element(genotype_posteriors.begin(), genotype_posteriors.end(), [] (const auto& lhs, const auto& rhs) { return lhs.second < rhs.second; });
        print_alleles(it->first);
        std::cout << std::endl;
        std::cout << it->second << std::endl;
        
//        for (const auto& scf : cf) {
//            auto n = scf.second[0] + scf.second[2] + scf.second[2];
//            std::cout << (scf.second[0] / n) << " " << (scf.second[1] / n) << " " << (scf.second[2] / n) << std::endl;
//        }
        
        return result;
    }
    
} // namespace Octopus
