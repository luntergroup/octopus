//
//  genotype_model.hpp
//  Octopus
//
//  Created by Daniel Cooke on 20/07/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_genotype_model_hpp
#define Octopus_genotype_model_hpp

#include <string>
#include <vector>
#include <unordered_map>

#include "common.hpp"
#include "mappable_map.hpp"
#include "aligned_read.hpp"
#include "haplotype.hpp"
#include "genotype.hpp"
#include "maths.hpp"

namespace Octopus
{
    namespace GenotypeModel
    {
        using ReadMap = MappableMap<SampleIdType, AlignedRead>;
        
        using HaplotypePriorCounts = std::unordered_map<Haplotype, double>;
        using HaplotypeFrequencies = std::unordered_map<Haplotype, double>;
        using GenotypeMarginals    = std::unordered_map<Genotype<Haplotype>, double>;
        using GenotypeLikelihoods  = std::unordered_map<SampleIdType, std::unordered_map<Genotype<Haplotype>, double>>;
        
        template <typename Genotype>
        double log_hardy_weinberg(const Genotype& genotype, const HaplotypeFrequencies& haplotype_frequencies)
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
        
        inline HaplotypeFrequencies init_haplotype_frequencies(const std::vector<Haplotype>& haplotypes)
        {
            HaplotypeFrequencies result {};
            result.reserve(haplotypes.size());
            
            const double uniform {1.0 / haplotypes.size()};
            
            for (const auto& haplotype : haplotypes) {
                result.emplace(haplotype, uniform);
            }
            
            return result;
        }
        
        inline double max_haplotype_frequency_change(const HaplotypeFrequencies& old_frequencies,
                                                     const HaplotypeFrequencies& new_frequencies)
        {
            double result {};
            
            for (const auto& h : new_frequencies) {
                auto change = std::abs(h.second - old_frequencies.at(h.first));
                if (change > result) result = change;
            }
            
            return result;
        }
        
    } // namespace GenotypeModel
} // namespace Octopus

#endif
