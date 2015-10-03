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
#include <algorithm> // std::transform
#include <iterator>  // std::begin, std::cbegin, std::cend, std::inserter
#include <cmath>     // std::abs

#include "common.hpp"
#include "mappable_map.hpp"
#include "aligned_read.hpp"
#include "haplotype.hpp"
#include "genotype.hpp"
#include "maths.hpp"
#include "reference_genome.hpp"
#include "haplotype_prior_model.hpp"

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
                auto num_occurences = genotype.count(haplotype);
                occurences.push_back(num_occurences);
                r += num_occurences * std::log(haplotype_frequencies.at(haplotype));
            }
            
            return Maths::log_multinomial_coefficient<double>(occurences) + r;
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
        
        inline HaplotypeFrequencies init_haplotype_frequencies(const HaplotypePriorCounts& haplotype_counts)
        {
            HaplotypeFrequencies result {};
            result.reserve(haplotype_counts.size());
            
            auto n = Maths::sum_values(haplotype_counts);
            
            for (const auto& haplotype_count : haplotype_counts) {
                result.emplace(haplotype_count.first, haplotype_count.second / n);
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
        
        inline HaplotypePriorCounts compute_haplotype_prior_counts(const std::vector<Haplotype>& haplotypes,
                                                                   ReferenceGenome& reference,
                                                                   HaplotypePriorModel& haplotype_prior_model)
        {
            using std::begin; using std::cbegin; using std::cend; using std::transform;
            
            HaplotypePriorCounts result {};
            
            if (haplotypes.empty()) return result;
            
            const Haplotype reference_haplotype {reference, haplotypes.front().get_region()};
            
            std::vector<double> p(haplotypes.size());
            transform(cbegin(haplotypes), cend(haplotypes), begin(p),
                      [&haplotype_prior_model, &reference_haplotype] (const auto& haplotype) {
                          return haplotype_prior_model.evaluate(haplotype, reference_haplotype);
                      });
            
            constexpr double precision {1'000.0};
            constexpr unsigned max_iterations {20};
            const auto alphas = Maths::maximum_likelihood_dirichlet_params(p, precision, max_iterations);
            
            result.reserve(haplotypes.size());
            
            transform(cbegin(haplotypes), cend(haplotypes), cbegin(alphas), std::inserter(result, begin(result)),
                      [] (const auto& haplotype, double a) { return std::make_pair(haplotype, a); });
            
            return result;
        }
        
    } // namespace GenotypeModel
} // namespace Octopus

#endif
