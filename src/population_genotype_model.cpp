//
//  population_genotype_model.cpp
//  Octopus
//
//  Created by Daniel Cooke on 26/08/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "population_genotype_model.h"

#include <cmath>

#include "read_model.h"
#include "maths.h"

#include <iostream> // TEST

namespace Octopus
{
    // public methods
    
    PopulationGenotypeModel::PopulationGenotypeModel(unsigned num_samples, unsigned sample_ploidy,
                                                     unsigned max_em_iterations, double em_epsilon)
    :
    max_em_iterations_ {max_em_iterations},
    em_epsilon_ {em_epsilon},
    num_samples_ {num_samples},
    sample_ploidy_ {sample_ploidy}
    {}
    
    // non member methods
    
    using HaplotypeFrequencies = std::unordered_map<Haplotype, double>;
    using GenotypeLikelihoods  = std::unordered_map<PopulationGenotypeModel::SampleIdType,
                                                std::unordered_map<Genotype<Haplotype>, double>>;
    using GenotypeMarginals    = std::unordered_map<Genotype<Haplotype>, double>;
    
    HaplotypeFrequencies init_haplotype_frequencies(const std::vector<Haplotype>& haplotypes)
    {
        HaplotypeFrequencies result {};
        result.reserve(haplotypes.size());
        
        double uniform {1.0 / haplotypes.size()};
        
        for (const auto& haplotype : haplotypes) {
            result.emplace(haplotype, uniform);
        }
        
        return result;
    }
    
    double log_hardy_weinberg(const Genotype<Haplotype>& genotype,
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
        
        auto c = log_multinomial_coefficient<double>(occurences.cbegin(), occurences.cend());
        
        return c * r;
    }
    
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
                                     const PopulationGenotypeModel::ReadMap& reads, unsigned ploidy)
    {
        ReadModel read_model {ploidy};
        
        GenotypeLikelihoods result {};
        result.reserve(reads.size());
        
        for (const auto& sample_reads : reads) {
            const auto& sample = sample_reads.first;
            std::unordered_map<Genotype<Haplotype>, double> likelihoods {};
            likelihoods.reserve(genotypes.size());
            
            for (const auto& genotype : genotypes) {
                likelihoods.emplace(genotype, read_model.log_probability(sample_reads.second.cbegin(),
                                                                        sample_reads.second.cend(), genotype, sample));
            }
            
            result.emplace(sample, std::move(likelihoods));
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
    
    template <typename K, typename V>
    std::vector<V> get_values(const std::unordered_map<K, V>& map)
    {
        std::vector<V> result(map.size());
        std::transform(std::cbegin(map), std::cend(map), std::begin(result),
                       [] (const auto& p) { return p.second; });
        return result;
    }
    
    void normalise(PopulationGenotypeModel::SampleGenotypeProbabilities& unnormalised_log_posteriors)
    {
        auto log_posteriors = get_values(unnormalised_log_posteriors);
        auto norm = log_sum_exp<double>(std::cbegin(log_posteriors), std::cend(log_posteriors));
        for (auto& p : unnormalised_log_posteriors) p.second -= norm;
    }
    
    void exponentiate(PopulationGenotypeModel::SampleGenotypeProbabilities& log_posteriors)
    {
        for (auto& p : log_posteriors) p.second = std::exp(p.second);
    }
    
    PopulationGenotypeModel::GenotypeProbabilities
    init_genotype_posteriors(const GenotypeMarginals& log_marginals,
                             const GenotypeLikelihoods& log_likilhoods)
    {
        PopulationGenotypeModel::GenotypeProbabilities result {};
        result.reserve(log_likilhoods.size());
        
        for (const auto& sample_likelihoods : log_likilhoods) {
            PopulationGenotypeModel::SampleGenotypeProbabilities gps {};
            gps.reserve(log_marginals.size());
            
            for (const auto& marginal : log_marginals) {
                gps.emplace(marginal.first, marginal.second + sample_likelihoods.second.at(marginal.first));
            }
            
            result.emplace(sample_likelihoods.first, std::move(gps));
        }
        
        for (auto& sample_genotypes : result) {
            normalise(sample_genotypes.second);
            exponentiate(sample_genotypes.second);
        }
        
        return result;
    }
    
    void update_posteriors(PopulationGenotypeModel::GenotypeProbabilities& genotype_posteriors,
                           const HaplotypeFrequencies& haplotype_frequencies,
                           GenotypeMarginals& marginal_genotype_log_probabilities,
                           const GenotypeLikelihoods& genotype_log_likilhoods)
    {
        update_marginal_genotype_log_probabilities(marginal_genotype_log_probabilities, haplotype_frequencies);
        
        for (auto& sample_genotype_log_likilhoods : genotype_posteriors) {
            auto sample = sample_genotype_log_likilhoods.first;
            for (auto& p : sample_genotype_log_likilhoods.second) {
                p.second = marginal_genotype_log_probabilities.at(p.first) + genotype_log_likilhoods.at(sample).at(p.first);
            }
        }
        
        for (auto& p : genotype_posteriors) {
            normalise(p.second);
            exponentiate(p.second);
        }
    }
    
    HaplotypeFrequencies
    compute_haplotype_frequencies(const PopulationGenotypeModel::GenotypeProbabilities& posteriors,
                                  unsigned num_samples, unsigned ploidy)
    {
        std::unordered_map<Haplotype, double> result {};
        
        for (const auto& sample_posteriors : posteriors) {
            for (const auto& posterior : sample_posteriors.second) {
                for (auto haplotype : posterior.first.get_unique()) {
                    result[haplotype] += posterior.first.num_occurences(haplotype) * posterior.second;
                }
            }
        }
        
        auto norm = num_samples * ploidy;
        
        for (auto& h : result) {
            h.second /= norm;
        }
        
        return result;
    }
    
    double max_haplotype_frequency_change(const HaplotypeFrequencies& old_frequencies,
                                          const HaplotypeFrequencies& new_frequencies)
    {
        double result {};
        
        for (const auto& h : new_frequencies) {
            auto change = std::abs(h.second - old_frequencies.at(h.first));
            if (change > result) result = change;
        }
        
        return result;
    }
    
    double do_em_iteration(PopulationGenotypeModel::GenotypeProbabilities& genotype_posteriors,
                           HaplotypeFrequencies& haplotype_frequencies,
                           GenotypeMarginals& marginal_genotype_log_probabilities,
                           const GenotypeLikelihoods& genotype_log_likilhoods,
                           unsigned ploidy)
    {
        update_posteriors(genotype_posteriors, haplotype_frequencies, marginal_genotype_log_probabilities, genotype_log_likilhoods);
        
        auto num_samples = static_cast<unsigned>(genotype_posteriors.size());
        
        auto new_frequencies = compute_haplotype_frequencies(genotype_posteriors, num_samples, ploidy);
        
        auto max_change = max_haplotype_frequency_change(haplotype_frequencies, new_frequencies);
        
        haplotype_frequencies = new_frequencies;
        
        return max_change;
    }
    
    // private methods
    
    PopulationGenotypeModel::GenotypeProbabilities
    PopulationGenotypeModel::do_evaluate(const std::vector<Haplotype>& haplotypes, const ReadMap& reads)
    {
        auto genotypes = get_all_genotypes(haplotypes, sample_ploidy_);
        
        std::cout << "there are " << genotypes.size() << " genotypes" << std::endl;
        
        const auto genotype_log_likilhoods = compute_genotype_log_likelihoods(genotypes, reads, sample_ploidy_);
        
        auto haplotype_frequencies               = init_haplotype_frequencies(haplotypes);
        auto marginal_genotype_log_probabilities = init_marginal_genotype_log_probabilities(genotypes, haplotype_frequencies);
        
        std::cout << "computed genotype log likilhoods" << std::endl;
        
        auto result = init_genotype_posteriors(marginal_genotype_log_probabilities, genotype_log_likilhoods);
        
        for (unsigned n {0}; n < max_em_iterations_; ++n) {
            std::cout << "EM iteration " << n << std::endl;
            auto c = do_em_iteration(result, haplotype_frequencies, marginal_genotype_log_probabilities,
                                     genotype_log_likilhoods, sample_ploidy_);
            if (c < em_epsilon_) break;
        }
        
//        for (const auto& h : haplotype_frequencies) {
//            std::cout << h.first << " " << h.second << std::endl;
//        }
        
        return result;
    }
    
} // end namespace Octopus

