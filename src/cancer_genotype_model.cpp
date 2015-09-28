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
    
    using HaplotypeFrequencies  = std::unordered_map<Haplotype, double>;
    using GenotypeMarginals     = std::unordered_map<Genotype<Haplotype>, double>;
    using ComponentPseudoCounts = std::unordered_map<Octopus::SampleIdType, std::array<double, 3>>;
    using ComponentFrequencies  = std::unordered_map<Octopus::SampleIdType, std::array<double, 3>>;
    
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
            r += num_occurences * haplotype_frequencies.at(haplotype);
        }
        
        return log_multinomial_coefficient<double>(occurences.cbegin(), occurences.cend()) * r;
    }
    
    double log_probability(const Genotype<Haplotype>& genotype, const std::array<double, 3>& f,
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
    
    double log_joint_liklihood(const Genotype<Haplotype>& g, const ComponentFrequencies& fs,
                               const ReadMap& reads, const HaplotypeFrequencies& pi,
                               const ComponentPseudoCounts& alphas)
    {
        double p {};
        
        for (const auto sample_reads : reads) {
            auto f = fs.at(sample_reads.first);
            auto a = alphas.at(sample_reads.first);
            p += log_dirichlet<double>(a.begin(), a.end(), f.begin()) + log_probability(g, f, sample_reads.second);
        }
        
        return log_hardy_weinberg(g, pi) + p;
    }
    
    GenotypeMarginals
    genotype_posteriors(const std::vector<Genotype<Haplotype>>& genotypes, const ReadMap& reads,
                        const HaplotypeFrequencies& pi, const ComponentPseudoCounts& alphas)
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
                p += log_probability(g, a, sample_reads.second);
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
    compute_haplotype_frequencies(const GenotypeMarginals& genotype_posteriors, unsigned ploidy)
    {
        HaplotypeFrequencies result {};
        
        for (const auto& posterior : genotype_posteriors) {
            for (auto haplotype : posterior.first.get_unique()) {
                result[haplotype] += posterior.first.num_occurences(haplotype) * posterior.second;
            }
        }
        
        for (auto& h : result) {
            h.second /= ploidy;
        }
        
        return result;
    }
    
    ComponentPseudoCounts
    compute_component_pseudo_counts(const ComponentPseudoCounts& prior_counts,
                                    const GenotypeMarginals& genotype_posteriors,
                                    const ReadMap& reads)
    {
        ComponentPseudoCounts result {};
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
    
    std::vector<std::array<double, 3>> generate_all_fs(unsigned n)
    {
        auto ratios = generate_all_ratios(n);
        std::vector<std::array<double, 3>> result(ratios.size());
        double epsilon {std::nextafter(0, 1.0)};
        std::transform(std::cbegin(ratios), std::cend(ratios), std::begin(result),
                       [n, epsilon] (const auto& arr) {
            return std::array<double, 3> {
                static_cast<double>(arr[0]) / n + epsilon,
                static_cast<double>(arr[1]) / n + epsilon,
                static_cast<double>(arr[2]) / n + epsilon
            };
        });
        return result;
    }
    
    // private methods
    
    CancerGenotypeModel::GenotypeProbabilities
    CancerGenotypeModel::do_evaluate(const std::vector<Haplotype>& haplotypes, const ReadMap& reads)
    {
        auto ratios = generate_all_fs(20);
        
//        for (auto r : ratios) {
//            std::cout << r[0] << " " << r[1] << " " << r[2] << std::endl;
//        }
        
        GenotypeProbabilities result {};
        
        ComponentPseudoCounts alphas {};
        for (const auto& s : reads) {
            alphas.emplace(s.first, std::array<double, 3> {1.0, 1.0, 1.0});
        }
        
        ComponentFrequencies fs {};
        for (const auto& s : reads) {
            fs.emplace(s.first, std::array<double, 3> {0.5, 0.5, std::nextafter(0, 1.0)});
        }
        
        HaplotypeFrequencies pi {};
        for (const auto& h : haplotypes) {
            pi.emplace(h, 1.0 / haplotypes.size());
        }
        
        auto genotypes = generate_all_genotypes(haplotypes, 3);
        
        remove_redundant_genotypes(genotypes);
        
        Genotype<Haplotype> gm {}; std::array<double, 3> fm {}; double m {-10000000};
        
        for (auto g : genotypes) {
            for (auto r : ratios) {
                ComponentFrequencies f {{reads.cbegin()->first, r}};
                auto h = log_joint_liklihood(g, f, reads, pi, alphas);
                if (h > m) {
                    gm = g;
                    fm = r;
                    m  = h;
                }
            }
        }
        
        print_alleles(gm);
        std::cout << std::endl;
        std::cout << fm[0] << " " << fm[1] << " " << fm[2] << std::endl;
        std::cout << m << std::endl;
        
        //std::cout << log_joint_liklihood(gm, fs, reads, pi, alphas) << std::endl;
        
        exit(0);
        
        //for (auto g : genotypes) std::cout << g << std::endl;
        
//        std::unordered_map<Genotype<Haplotype>, double> ll(genotypes.size());
//        for (const auto& genotype : genotypes) {
//            ll.emplace(genotype, log_joint_liklihood(genotype, fs, reads, pi, alphas));
//        }
        
        auto gp = genotype_posteriors(genotypes, reads, pi, alphas);
        
        auto hf = compute_haplotype_frequencies(gp, 3);
        auto cf = compute_component_pseudo_counts(alphas, gp, reads);
        
//        for (const auto f : hf) {
//            std::cout << f.second << " ";
//        }
//        std::cout << std::endl;
        
//        for (const auto g : ll) {
//            std::cout << g << std::endl;
//        }
        
        auto it = std::max_element(gp.begin(), gp.end(), [] (const auto& lhs, const auto& rhs) { return lhs.second < rhs.second; });
        //std::cout << it->first << std::endl;
        print_alleles(it->first);
        std::cout << std::endl;
        std::cout << it->second << std::endl;
        
        for (const auto& scf : cf) {
            auto n = scf.second[0] + scf.second[2] + scf.second[2];
            std::cout << (scf.second[0] / n) << " " << (scf.second[1] / n) << " " << (scf.second[2] / n) << std::endl;
        }
        
        return result;
    }
    
} // end namespace Octopus
