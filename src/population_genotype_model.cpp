//
//  population_genotype_model.cpp
//  Octopus
//
//  Created by Daniel Cooke on 26/08/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "population_genotype_model.hpp"

#include <cmath>

#include "read_model.hpp"
#include "maths.hpp"

#include <iostream> // TEST

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
    
    GenotypeMarginals
    init_marginal_genotype_log_probabilities(const std::vector<Genotype<Haplotype>>& genotypes,
                                             const HaplotypeFrequencies& haplotype_frequencies)
    {
        GenotypeMarginals result {};
        result.reserve(genotypes.size());
        
        for (const auto& genotype : genotypes) {
            result.emplace(genotype, log_hardy_weinberg(genotype, haplotype_frequencies));
        }
        
        return result;
    }
    
    GenotypeLikelihoods
    compute_genotype_log_likelihoods(const std::vector<Genotype<Haplotype>>& genotypes,
                                     const ReadMap& reads)
    {
        ReadModel read_model {genotypes.front().ploidy()};
        
        GenotypeLikelihoods result {};
        result.reserve(reads.size());
        
        for (const auto& sample_reads : reads) {
            const auto& sample = sample_reads.first;
            std::unordered_map<Genotype<Haplotype>, double> likelihoods {};
            likelihoods.reserve(genotypes.size());
            
            for (const auto& genotype : genotypes) {
                likelihoods.emplace(genotype, read_model.log_probability(std::cbegin(sample_reads.second),
                                                                         std::cend(sample_reads.second),
                                                                         genotype));
            }
            
            result.emplace(sample, std::move(likelihoods));
        }
        
        return result;
    }
    
    template <typename K, typename V>
    static std::vector<V> get_values(const std::unordered_map<K, V>& map)
    {
        std::vector<V> result(map.size());
        std::transform(std::cbegin(map), std::cend(map), std::begin(result),
                       [] (const auto& p) { return p.second; });
        return result;
    }
    
    static void normalise(Population::SampleGenotypeProbabilities& unnormalised_log_posteriors)
    {
        auto log_posteriors = get_values(unnormalised_log_posteriors);
        const auto norm = Maths::log_sum_exp<double>(log_posteriors);
        for (auto& p : unnormalised_log_posteriors) p.second -= norm;
    }
    
    static void exponentiate(Population::SampleGenotypeProbabilities& log_posteriors)
    {
        for (auto& p : log_posteriors) p.second = std::exp(p.second);
    }
    
    Population::GenotypeProbabilities
    init_genotype_posteriors(const GenotypeMarginals& log_marginals,
                             const GenotypeLikelihoods& log_likilhoods)
    {
        Population::GenotypeProbabilities result {};
        result.reserve(log_likilhoods.size());
        
        for (const auto& sample_likelihoods : log_likilhoods) {
            Population::SampleGenotypeProbabilities gps {};
            gps.reserve(log_marginals.size());
            
            for (const auto& log_marginal : log_marginals) { // (genotype, log prob)
                gps.emplace(log_marginal.first, log_marginal.second + sample_likelihoods.second.at(log_marginal.first));
            }
            
            normalise(gps);
            exponentiate(gps);
            
            result.emplace(sample_likelihoods.first, std::move(gps));
        }
        
        return result;
    }
    
    void update_marginal_genotype_log_probabilities(GenotypeMarginals& marginals,
                                                    const HaplotypeFrequencies& haplotype_frequencies)
    {
        for (auto& p : marginals) {
            p.second = log_hardy_weinberg(p.first, haplotype_frequencies);
        }
    }
    
    void update_genotype_posteriors(Population::GenotypeProbabilities& genotype_posteriors,
                                    const HaplotypeFrequencies& haplotype_frequencies,
                                    const GenotypeMarginals& marginal_genotype_log_probabilities,
                                    const GenotypeLikelihoods& genotype_log_likilhoods)
    {
        for (auto& sample_genotype_posteriors : genotype_posteriors) {
            const auto& sample = sample_genotype_posteriors.first;
            
            for (auto& gp : sample_genotype_posteriors.second) { // (genotype, probability)
                gp.second = marginal_genotype_log_probabilities.at(gp.first) + genotype_log_likilhoods.at(sample).at(gp.first);
            }
            
            normalise(sample_genotype_posteriors.second);
            exponentiate(sample_genotype_posteriors.second);
        }
    }
    
    HaplotypeFrequencies
    compute_haplotype_frequencies(const std::vector<Haplotype>& haplotypes,
                                  const std::vector<Genotype<Haplotype>>& genotypes,
                                  const Population::GenotypeProbabilities& genotype_posteriors)
    {
        HaplotypeFrequencies result {};
        result.reserve(haplotypes.size());
        
        for (const auto& haplotype : haplotypes) {
            double p {};
            
            for (const auto& genotype : genotypes) {
                auto n = genotype.count(haplotype);
                if (n > 0) {
                    double t {};
                    for (const auto& sample_genotype_posteriors : genotype_posteriors) {
                        t += sample_genotype_posteriors.second.at(genotype);
                    }
                    p += n * t;
                }
            }
            
            result.emplace(haplotype, p);
        }
        
        const auto norm = genotype_posteriors.size() * genotypes.front().ploidy();
        
        for (auto& hf : result) {
            hf.second /= norm;
        }
        
        return result;
    }
    
    HaplotypeFrequencies
    compute_haplotype_frequencies(const HaplotypePriorCounts& haplotype_prior_counts,
                                  const std::vector<Genotype<Haplotype>>& genotypes,
                                  const Population::GenotypeProbabilities& genotype_posteriors)
    {
        HaplotypeFrequencies result {};
        result.reserve(haplotype_prior_counts.size());
        
        for (const auto& haplotype_count : haplotype_prior_counts) {
            auto p = haplotype_count.second;
            
            for (const auto& genotype : genotypes) {
                auto n = genotype.count(haplotype_count.first);
                if (n > 0) {
                    double t {};
                    for (const auto& sample_genotype_posteriors : genotype_posteriors) {
                        t += sample_genotype_posteriors.second.at(genotype);
                    }
                    p += n * t;
                }
            }
            
            result.emplace(haplotype_count.first, p);
        }
        
        const auto norm = genotype_posteriors.size() * genotypes.front().ploidy() + Maths::sum_values(haplotype_prior_counts);
        
        for (auto& hf : result) {
            hf.second /= norm;
        }
        
        return result;
    }
    
    double do_em_iteration(const std::vector<Haplotype>& haplotypes,
                           const std::vector<Genotype<Haplotype>>& genotypes,
                           Population::GenotypeProbabilities& genotype_posteriors,
                           HaplotypeFrequencies& haplotype_frequencies,
                           GenotypeMarginals& marginal_genotype_log_probabilities,
                           const GenotypeLikelihoods& genotype_log_likilhoods,
                           const HaplotypePriorCounts& haplotype_prior_counts)
    {
        auto new_frequencies = compute_haplotype_frequencies(haplotype_prior_counts, genotypes, genotype_posteriors);
        
        auto max_change = max_haplotype_frequency_change(haplotype_frequencies, new_frequencies);
        
        update_marginal_genotype_log_probabilities(marginal_genotype_log_probabilities, new_frequencies);
        update_genotype_posteriors(genotype_posteriors, new_frequencies,
                                   marginal_genotype_log_probabilities, genotype_log_likilhoods);
        
        haplotype_frequencies = new_frequencies;
        
        return max_change;
    }
    
    // private methods
    
        void print_top_haplotypes(const std::vector<Haplotype>& haplotypes, const ReadMap& reads, size_t n = 3)
        {
            std::cout << "top " << n << " haplotype likelihoods for each read in each sample" << std::endl;
            auto m = std::min(n, haplotypes.size());
            for (const auto& sample_reads : reads) {
                std::cout << sample_reads.first << ":" << std::endl;
                for (const auto& read : sample_reads.second) {
                    std::cout << read.get_region() << " " << read.get_cigar_string() << ":" << std::endl;
                    std::vector<std::pair<Haplotype, double>> top {};
                    top.reserve(haplotypes.size());
                    for (const auto& haplotype : haplotypes) {
                        top.emplace_back(haplotype, ReadModel(2).log_probability(read, haplotype));
                    }
                    std::sort(std::begin(top), std::end(top), [] (const auto& lhs, const auto& rhs) {
                        return lhs.second > rhs.second;
                    });
                    for (unsigned i {}; i < m; ++i) {
                        std::cout << "\t* ";
                        print_variant_alleles(top[i].first);
                        std::cout << " " << top[i].second << std::endl;
                    }
                }
            }
        }
        
        void print_top_haplotypes(const std::vector<Haplotype>& haplotypes, const ReadMap& reads,
                                  const SampleIdType& sample, size_t n = 3)
        {
            std::cout << "top " << n << " haplotype likelihoods for each read in sample " << sample << std::endl;
            auto m = std::min(n, haplotypes.size());
            for (const auto& read : reads.at(sample)) {
                std::cout << read.get_region() << " " << read.get_cigar_string() << ":" << std::endl;
                std::vector<std::pair<Haplotype, double>> top {};
                top.reserve(haplotypes.size());
                for (const auto& haplotype : haplotypes) {
                    top.emplace_back(haplotype, ReadModel(2).log_probability(read, haplotype));
                }
                std::sort(std::begin(top), std::end(top), [] (const auto& lhs, const auto& rhs) {
                    return lhs.second > rhs.second;
                });
                for (unsigned i {}; i < m; ++i) {
                    std::cout << "\t* ";
                    print_variant_alleles(top[i].first);
                    std::cout << " " << top[i].second << std::endl;
                }
            }
        }
        
        void print_top_genotypes(const std::vector<Genotype<Haplotype>>& genotypes, const ReadMap& reads, size_t n = 3)
        {
            std::cout << "top " << n << " genotype likelihoods for each read in each sample" << std::endl;
            auto m = std::min(n, genotypes.size());
            for (const auto& sample_reads : reads) {
                std::cout << sample_reads.first << ":" << std::endl;
                for (const auto& read : sample_reads.second) {
                    std::cout << read.get_region() << " " << read.get_cigar_string() << ":" << std::endl;
                    std::vector<std::pair<Genotype<Haplotype>, double>> top {};
                    top.reserve(genotypes.size());
                    for (const auto& genotype : genotypes) {
                        top.emplace_back(genotype, ReadModel(2).log_probability(read, genotype));
                    }
                    std::sort(std::begin(top), std::end(top), [] (const auto& lhs, const auto& rhs) {
                        return lhs.second > rhs.second;
                    });
                    for (unsigned i {}; i < m; ++i) {
                        std::cout << "\t* ";
                        print_variant_alleles(top[i].first);
                        std::cout << " " << top[i].second << std::endl;
                    }
                }
            }
        }
        
        void print_top_genotypes(const std::vector<Genotype<Haplotype>>& genotypes, const ReadMap& reads,
                                 const SampleIdType& sample, size_t n = 3)
        {
            std::cout << "top " << n << " genotype likelihoods for each read in sample " << sample << std::endl;
            auto m = std::min(n, genotypes.size());
            for (const auto& read : reads.at(sample)) {
                std::cout << read.get_region() << " " << read.get_cigar_string() << ":" << std::endl;
                std::vector<std::pair<Genotype<Haplotype>, double>> top {};
                top.reserve(genotypes.size());
                for (const auto& genotype : genotypes) {
                    top.emplace_back(genotype, ReadModel(2).log_probability(read, genotype));
                }
                std::sort(std::begin(top), std::end(top), [] (const auto& lhs, const auto& rhs) {
                    return lhs.second > rhs.second;
                });
                for (unsigned i {}; i < m; ++i) {
                    std::cout << "\t* ";
                    print_variant_alleles(top[i].first);
                    std::cout << " " << top[i].second << std::endl;
                }
            }
        }
        
    Population::Latents
    Population::evaluate(const std::vector<Haplotype>& haplotypes, const ReadMap& reads, ReferenceGenome& reference)
    {
        auto genotypes = generate_all_genotypes(haplotypes, ploidy_);
        
        const auto genotype_log_likilhoods = compute_genotype_log_likelihoods(genotypes, reads);
        
//        print_top_haplotypes(haplotypes, reads, "HG00102", 10);
//        std::cout << std::endl;
//        print_top_genotypes(genotypes, reads, "HG00102", 10);
//        std::cout << std::endl;
//        print(genotype_log_likilhoods, 10);
//        std::cout << std::endl;
        //exit(0);
        
        auto haplotype_prior_counts = compute_haplotype_prior_counts(haplotypes, reference, haplotype_prior_model_);
        
        auto haplotype_frequencies               = init_haplotype_frequencies(haplotype_prior_counts);
        auto marginal_genotype_log_probabilities = init_marginal_genotype_log_probabilities(genotypes, haplotype_frequencies);
        
        auto genotype_posteriors = init_genotype_posteriors(marginal_genotype_log_probabilities, genotype_log_likilhoods);
        
//        print(haplotype_prior_counts);
//        std::cout << std::endl;
//        print(genotype_posteriors, 10);
//        std::cout << std::endl;
        //exit(0);
        
        for (unsigned n {0}; n < max_em_iterations_; ++n) {
            // std::cout << "EM iteration " << n << std::endl;
            auto c = do_em_iteration(haplotypes, genotypes, genotype_posteriors, haplotype_frequencies,
                                     marginal_genotype_log_probabilities, genotype_log_likilhoods,
                                     haplotype_prior_counts);
            if (c < em_epsilon_) break;
        }
        
        //print(genotype_posteriors, 10);
        //exit(0);
        
        return Latents {std::move(genotype_posteriors), std::move(haplotype_frequencies)};
    }
    
    } // namespace GenotypeModel
} // namespace Octopus

