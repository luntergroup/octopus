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
                likelihoods.emplace(genotype, read_model.log_probability(sample_reads.second.cbegin(),
                                                                         sample_reads.second.cend(),
                                                                         genotype, sample));
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
    
    void normalise(Population::SampleGenotypeProbabilities& unnormalised_log_posteriors)
    {
        auto log_posteriors = get_values(unnormalised_log_posteriors);
        const auto norm = log_sum_exp<double>(std::cbegin(log_posteriors), std::cend(log_posteriors));
        for (auto& p : unnormalised_log_posteriors) p.second -= norm;
    }
    
    void exponentiate(Population::SampleGenotypeProbabilities& log_posteriors)
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
            
            for (const auto& marginal : log_marginals) {
                gps.emplace(marginal.first, marginal.second + sample_likelihoods.second.at(marginal.first));
            }
            
            normalise(gps);
            exponentiate(gps);
            
            result.emplace(sample_likelihoods.first, std::move(gps));
        }
        
        return result;
    }
    
    void update_genotype_posteriors(Population::GenotypeProbabilities& genotype_posteriors,
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
    compute_haplotype_frequencies(const std::vector<Haplotype>& haplotypes,
                                  const std::vector<Genotype<Haplotype>>& genotypes,
                                  const Population::GenotypeProbabilities& genotype_posteriors)
    {
        HaplotypeFrequencies result {};
        result.reserve(haplotypes.size());
        
        for (const auto& haplotype : haplotypes) {
            double p {};
            
            for (const auto& genotype : genotypes) {
                auto n = genotype.num_occurences(haplotype);
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
    
    double do_em_iteration(const std::vector<Haplotype>& haplotypes,
                           const std::vector<Genotype<Haplotype>>& genotypes,
                           Population::GenotypeProbabilities& genotype_posteriors,
                           HaplotypeFrequencies& haplotype_frequencies,
                           GenotypeMarginals& marginal_genotype_log_probabilities,
                           const GenotypeLikelihoods& genotype_log_likilhoods)
    {
        auto new_frequencies = compute_haplotype_frequencies(haplotypes, genotypes, genotype_posteriors);
        
        update_genotype_posteriors(genotype_posteriors, haplotype_frequencies,
                                   marginal_genotype_log_probabilities, genotype_log_likilhoods);
        
        auto max_change = max_haplotype_frequency_change(haplotype_frequencies, new_frequencies);
        
        update_genotype_posteriors(genotype_posteriors, new_frequencies,
                                   marginal_genotype_log_probabilities, genotype_log_likilhoods);
        
        haplotype_frequencies = new_frequencies;
        
        return max_change;
    }
    
    // private methods
    
    Population::Latents Population::evaluate(const std::vector<Haplotype>& haplotypes, const ReadMap& reads)
    {
        auto genotypes = generate_all_genotypes(haplotypes, ploidy_);
        
        //std::cout << "there are " << genotypes.size() << " genotypes" << std::endl;
        
        auto haplotype_frequencies               = init_haplotype_frequencies(haplotypes);
        auto marginal_genotype_log_probabilities = init_marginal_genotype_log_probabilities(genotypes, haplotype_frequencies);
        
        const auto genotype_log_likilhoods = compute_genotype_log_likelihoods(genotypes, reads);
        
        auto genotype_posteriors = init_genotype_posteriors(marginal_genotype_log_probabilities, genotype_log_likilhoods);
        
        for (unsigned n {0}; n < max_em_iterations_; ++n) {
            //std::cout << "EM iteration " << n << std::endl;
            auto c = do_em_iteration(haplotypes, genotypes, genotype_posteriors, haplotype_frequencies,
                                     marginal_genotype_log_probabilities, genotype_log_likilhoods);
            if (c < em_epsilon_) break;
        }
        
        Latents result {};
        result.genotype_posteriors   = std::move(genotype_posteriors);
        result.haplotype_frequencies = std::move(haplotype_frequencies);
        
        return result;
    }
    
    } // namespace GenotypeModel
} // namespace Octopus

