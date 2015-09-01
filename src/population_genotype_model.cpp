//
//  population_genotype_model.cpp
//  Octopus
//
//  Created by Daniel Cooke on 26/08/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "population_genotype_model.h"

#include "read_model.h"
#include "maths.h"

namespace Octopus
{
    // public methods
    
    PopulationGenotypeModel::PopulationGenotypeModel(unsigned num_samples, unsigned sample_ploidy)
    :
    num_samples_ {num_samples},
    sample_ploidy_ {sample_ploidy}
    {}
    
    // non member methods
    
    std::unordered_map<Haplotype, double> init_haplotype_frequencies(const std::vector<Haplotype>& haplotypes)
    {
        std::unordered_map<Haplotype, double> result {};
        result.reserve(haplotypes.size());
        
        double uniform {1.0 / haplotypes.size()};
        
        for (const auto& haplotype : haplotypes) {
            result.emplace(haplotype, uniform);
        }
        
        return result;
    }
    
    double log_hardy_weinberg(const Genotype<Haplotype>& genotype,
                              const std::unordered_map<Haplotype, double>& haplotype_frequencies)
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
    
    std::unordered_map<Genotype<Haplotype>, double>
    init_marginal_genotype_log_probabilities(const std::vector<Genotype<Haplotype>>& genotypes,
                                             const std::unordered_map<Haplotype, double>& haplotype_frequencies)
    {
        std::unordered_map<Genotype<Haplotype>, double> result {};
        result.reserve(genotypes.size());
        
        for (const auto& genotype : genotypes) {
            result.emplace(genotype, log_hardy_weinberg(genotype, haplotype_frequencies));
        }
        
        return result;
    }
    
    void update_marginal_genotype_log_probabilities(std::unordered_map<Genotype<Haplotype>, double>& marginals,
                                                    const std::unordered_map<Haplotype, double>& haplotype_frequencies)
    {
        for (auto& p : marginals) {
            p.second = log_hardy_weinberg(p.first, haplotype_frequencies);
        }
    }
    
    // private methods
    
    PopulationGenotypeModel::GenotypeProbabilities
    PopulationGenotypeModel::do_evaluate(const std::vector<Haplotype>& haplotypes, const ReadMap& reads)
    {
        auto genotypes = get_all_genotypes(haplotypes, sample_ploidy_);
        
        ReadModel read_model {sample_ploidy_};
        
        auto haplotype_frequencies = init_haplotype_frequencies(haplotypes);
        
        auto marginal_genotype_log_probabilities = init_marginal_genotype_log_probabilities(genotypes, haplotype_frequencies);
        
        GenotypeProbabilities result {};
        result.reserve(num_samples_);
        
        for (const auto& sample_reads : reads) {
            auto sample = sample_reads.first;
            
            SampleGenotypeProbabilities gps {};
            gps.reserve(genotypes.size());
            
            for (const auto& genotype : genotypes) {
                auto p = read_model.log_probability(sample_reads.second.cbegin(), sample_reads.second.cend(), genotype, sample);
                gps.emplace(genotype, marginal_genotype_log_probabilities[genotype] + p);
            }
            
            result.emplace(sample, std::move(gps));
        }
        
        // exp and normalise
        
        std::vector<double> log_probs(genotypes.size());
        for (auto& sample_genotypes : result) {
            std::transform(sample_genotypes.second.cbegin(), sample_genotypes.second.cend(), log_probs.begin(), [] (const auto& p) { return p.second; });
            auto n = log_sum_exp<double>(log_probs.cbegin(), log_probs.cend());
            for (auto& p : sample_genotypes.second) {
                p.second = std::exp(p.second - n);
            }
        }
        
        return result;
    }
    
} // end namespace Octopus

