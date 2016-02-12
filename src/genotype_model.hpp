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
#include <algorithm>
#include <numeric>
#include <iterator>
#include <cmath>

#include <iostream> // TEST

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
        using HaplotypePriorCountMap = std::unordered_map<Haplotype, double>;
        using HaplotypeFrequencyMap  = std::unordered_map<Haplotype, double>;
        
        namespace detail {
            template <typename Genotype>
            double log_hardy_weinberg_haploid(const Genotype& genotype,
                                              const HaplotypeFrequencyMap& haplotype_frequencies)
            {
                return std::log(haplotype_frequencies.at(genotype[0]));
            }
            
            template <typename Genotype>
            double log_hardy_weinberg_diploid(const Genotype& genotype,
                                              const HaplotypeFrequencyMap& haplotype_frequencies)
            {
                static const double ln_2 {std::log(2.0)};
                
                if (genotype.is_homozygous()) {
                    return 2 * std::log(haplotype_frequencies.at(genotype[0]));
                } else {
                    return std::log(haplotype_frequencies.at(genotype[0])) + std::log(haplotype_frequencies.at(genotype[1])) + ln_2;
                }
            }
            
            template <typename Genotype>
            double log_hardy_weinberg_triploid(const Genotype& genotype,
                                               const HaplotypeFrequencyMap& haplotype_frequencies)
            {
//                if (genotype.is_homozygous()) {
//                    return 3 * std::log(haplotype_frequencies.at(genotype[0]));
//                }
//                
                // TODO: other cases
                
                auto unique_haplotypes = genotype.copy_unique();
                
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
            
            template <typename Genotype>
            double log_hardy_weinberg_polyploid(const Genotype& genotype,
                                                const HaplotypeFrequencyMap& haplotype_frequencies)
            {
                auto unique_haplotypes = genotype.copy_unique();
                
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
        }
        
        // TODO: improve this, possible bottleneck in EM update at the moment
        template <typename Genotype>
        double log_hardy_weinberg(const Genotype& genotype, const HaplotypeFrequencyMap& haplotype_frequencies)
        {
            switch (genotype.ploidy()) {
                case 1 : return detail::log_hardy_weinberg_haploid(genotype, haplotype_frequencies);
                case 2 : return detail::log_hardy_weinberg_diploid(genotype, haplotype_frequencies);
                case 3 : return detail::log_hardy_weinberg_triploid(genotype, haplotype_frequencies);
                default: return detail::log_hardy_weinberg_polyploid(genotype, haplotype_frequencies);
            }
        }
        
        inline HaplotypeFrequencyMap init_haplotype_frequencies(const std::vector<Haplotype>& haplotypes)
        {
            HaplotypeFrequencyMap result {};
            result.reserve(haplotypes.size());
            
            const double uniform {1.0 / haplotypes.size()};
            
            for (const auto& haplotype : haplotypes) {
                result.emplace(haplotype, uniform);
            }
            
            return result;
        }
        
        inline HaplotypeFrequencyMap init_haplotype_frequencies(const HaplotypePriorCountMap& haplotype_counts)
        {
            HaplotypeFrequencyMap result {};
            result.reserve(haplotype_counts.size());
            
            auto n = Maths::sum_values(haplotype_counts);
            
            for (const auto& haplotype_count : haplotype_counts) {
                result.emplace(haplotype_count.first, haplotype_count.second / n);
            }
            
            return result;
        }
        
        inline HaplotypePriorCountMap compute_haplotype_prior_counts(const std::vector<Haplotype>& haplotypes,
                                                                   const ReferenceGenome& reference,
                                                                   HaplotypePriorModel& haplotype_prior_model)
        {
            using std::begin; using std::cbegin; using std::cend; using std::transform;
            
            HaplotypePriorCountMap result {};
            
            if (haplotypes.empty()) return result;
            
            const Haplotype reference_haplotype {haplotypes.front().get_region(), reference};
            
            std::vector<double> p(haplotypes.size());
            transform(cbegin(haplotypes), cend(haplotypes), begin(p),
                      [&haplotype_prior_model, &reference_haplotype] (const auto& haplotype) {
                          return haplotype_prior_model.evaluate(haplotype, reference_haplotype);
                      });
            
            const auto norm = std::accumulate(cbegin(p), cend(p), 0.0);
            
            transform(cbegin(p), cend(p), begin(p), [norm] (auto x) { return x / norm; });
            
            constexpr double precision {40.0};
            constexpr unsigned max_iterations {100};
            const auto alphas = Maths::dirichlet_mle(p, precision, max_iterations);
            
            result.reserve(haplotypes.size());
            
            transform(cbegin(haplotypes), cend(haplotypes), cbegin(alphas), std::inserter(result, begin(result)),
                      [] (const auto& haplotype, double a) { return std::make_pair(haplotype, a); });
            
            for (auto& h : result) h.second = 1; // DEBUG - uniform priors
            
            return result;
        }
    } // namespace GenotypeModel
} // namespace Octopus

#endif
