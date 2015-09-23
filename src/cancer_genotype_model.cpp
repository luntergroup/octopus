//
//  cancer_genotype_model.cpp
//  Octopus
//
//  Created by Daniel Cooke on 26/08/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "cancer_genotype_model.h"

#include <array>
#include <numeric>
#include <algorithm>

#include "read_model.h"
#include "common.h"
#include "maths.h"

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
    
    using FMap = std::unordered_map<Octopus::SampleIdType, std::array<double, 3>>;
    using HaplotypeFrequencies = std::unordered_map<Haplotype, double>;
    using GenotypeMarginals    = std::unordered_map<Genotype<Haplotype>, double>;
    
    double log_hw(const Genotype<Haplotype>& genotype, const HaplotypeFrequencies& haplotype_frequencies)
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
    
    double log_joint_liklihood(const Genotype<Haplotype>& g, const FMap& fs, const ReadMap& reads,
                               const HaplotypeFrequencies& pi, const FMap& alphas)
    {
        ReadModel rm {3};
        
        double p {};
        
        for (const auto sample_reads : reads) {
            auto f = fs.at(sample_reads.first);
            auto a = alphas.at(sample_reads.first);
            
            p += log_dirichlet<double>(a.begin(), a.end(), f.begin());
            
            for (const auto read : sample_reads.second) {
                p += log_sum_exp(std::log(f[0]) + rm.log_probability(read, g[0], sample_reads.first),
                                 std::log(f[1]) + rm.log_probability(read, g[1], sample_reads.first),
                                 std::log(f[2]) + rm.log_probability(read, g[2], sample_reads.first));
            }
        }
        
        return log_hw(g, pi) + p;
    }
    
    GenotypeMarginals
    genotype_posteriors(const std::vector<Genotype<Haplotype>>& genotypes, const ReadMap& reads,
                        const HaplotypeFrequencies& pi, const FMap& alphas)
    {
        ReadModel rm {3};
        
        GenotypeMarginals result {};
        result.reserve(genotypes.size());
        
        std::vector<double> t {};
        t.reserve(genotypes.size());
        
        for (const auto& g : genotypes) {
            auto p = log_hw(g, pi);
            
            for (const auto& sample_reads : reads) {
                const auto& s = sample_reads.first;
                const auto a = alphas.at(s);
                for (const auto& read : sample_reads.second) {
                    p += log_sum_exp(std::log(a[0]) + rm.log_probability(read, g[0], sample_reads.first),
                                     std::log(a[1]) + rm.log_probability(read, g[1], sample_reads.first),
                                     std::log(a[2]) + rm.log_probability(read, g[2], sample_reads.first));
                }
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
    
//    std::unordered_map<Octopus::SampleIdType, double>
//    fmap_posterior(const FMap& f, const std::vector<Genotype<Haplotype>>& genotypes, const ReadMap& reads,
//                   const HaplotypeFrequencies& pi, const FMap& alphas)
//    {
//        
//    }
    
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
    
    void remove_redundant_genotypes(std::vector<Genotype<Haplotype>>& genotypes)
    {
        genotypes.erase(std::remove_if(std::begin(genotypes), std::end(genotypes),
                                       [] (const auto& genotype) {
                                           return genotype.zygosity() != 3;
                                       }), std::end(genotypes));
    }
    
    // private methods
    
    CancerGenotypeModel::GenotypeProbabilities
    CancerGenotypeModel::do_evaluate(const std::vector<Haplotype>& haplotypes, const ReadMap& reads)
    {
        GenotypeProbabilities result {};
        
        FMap alphas {};
        for (const auto& s : reads) {
            alphas.emplace(s.first, std::array<double, 3> {1.0, 1.0, 1.0});
        }
        
        FMap fs {};
        for (const auto& s : reads) {
            fs.emplace(s.first, std::array<double, 3> {0.1, 0.1, 0.8});
        }
        
        HaplotypeFrequencies pi {};
        for (const auto& h : haplotypes) {
            pi.emplace(h, 1.0 / haplotypes.size());
        }
        
        auto genotypes = generate_all_genotypes(haplotypes, 3);
        
        //remove_redundant_genotypes(genotypes);
        
        //for (auto g : genotypes) std::cout << g << std::endl;
        
//        std::unordered_map<Genotype<Haplotype>, double> ll(genotypes.size());
//        for (const auto& genotype : genotypes) {
//            ll.emplace(genotype, log_joint_liklihood(genotype, fs, reads, pi, alphas));
//        }
        
        auto ll = genotype_posteriors(genotypes, reads, pi, alphas);
        
//        for (const auto g : ll) {
//            std::cout << g << std::endl;
//        }
        
        auto it = std::max_element(ll.begin(), ll.end(), [] (const auto& lhs, const auto& rhs) { return lhs.second < rhs.second; });
        std::cout << it->first << std::endl;
        std::cout << it->second << std::endl;
        
        return result;
    }
    
} // end namespace Octopus
