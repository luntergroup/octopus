//
//  bayesian_genotype_model.h
//  Octopus
//
//  Created by Daniel Cooke on 12/05/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_bayesian_genotype_model_h
#define Octopus_bayesian_genotype_model_h

#include <vector>
#include <string>
#include <unordered_map>
#include <algorithm>  // std::transform, std::copy_if, std::any_of, std::binary_search
#include <functional> // std::reference_wrapper
#include <cmath>      // std::exp, std::log
#include <cstddef>    // std::size_t
#include <numeric>    // std::accumulate

#include "haplotype.h"
#include "genotype.h"
#include "read_model.h"
#include "maths.h"
#include "aligned_read.h"
#include "pair_hmm.h"

namespace Octopus
{

namespace BayesianGenotypeModel
{
    template <typename RealType>
    using HaplotypePseudoCounts = std::unordered_map<Haplotype, RealType>;
    
    template <typename RealType>
    using SampleGenotypeProbabilities = std::unordered_map<Genotype<Haplotype>, RealType>;
    
    template <typename SampleIdType, typename RealType>
    using GenotypeProbabilities = std::unordered_map<SampleIdType, SampleGenotypeProbabilities<RealType>>;
    
    template <typename SampleIdType, typename RealType>
    struct Latents
    {
        HaplotypePseudoCounts<RealType> haplotype_pseudo_counts;
        GenotypeProbabilities<SampleIdType, RealType> genotype_posteriors;
        
        Latents() = default;
        template <typename HaplotypePseudoCounts_, typename GenotypeProbabilities_>
        Latents(HaplotypePseudoCounts_&& haplotype_pseudo_counts, GenotypeProbabilities_&& genotype_posteriors)
        :
        haplotype_pseudo_counts {std::forward<HaplotypePseudoCounts_>(haplotype_pseudo_counts)},
        genotype_posteriors {std::forward<GenotypeProbabilities_>(genotype_posteriors)}
        {}
    };
    
    template <typename MapType, typename RealType>
    HaplotypePseudoCounts<RealType>
    get_haplotype_prior_pseudo_counts(const MapType& the_haplotype_priors,
                                      const Haplotype& the_reference_haplotype,
                                      RealType the_reference_haplotype_pseudo_count=1)
    {
        HaplotypePseudoCounts<RealType> result {};
        result.reserve(the_haplotype_priors.size());
        
        RealType reference_bias {0.9};
        static const RealType concentration {std::log(7.0)};
        
        for (const auto& haplotype_prior : the_haplotype_priors) {
            result.emplace(haplotype_prior.first, digamma_inv(reference_bias * std::log(haplotype_prior.second) - concentration));
        }
        
//        std::cout << "prior pseudo counts: " << std::endl;
//        for (const auto& count : result) {
//            std::cout << count.second << std::endl;
//        }
//        std::cout << "sum: " << sum_values(result) << std::endl;
        
        result[the_reference_haplotype] += reference_bias;
        
        //for (auto& count : result) count.second = 0.5; // TEST
        
        return result;
    }
    
    template <typename RealType, typename Genotypes>
    RealType probability_haplotype_in_samples(const Haplotype& haplotype,
                                              const Genotypes& genotypes,
                                              const HaplotypePseudoCounts<RealType>& haplotype_pseudo_counts)
    {
        RealType result {0};
        
        // TODO: this calculation is incorrect
        for (const auto& genotype : genotypes) {
            if (genotype.contains(haplotype)) {
                result += posterior_predictive_probability(genotype, haplotype_pseudo_counts);
            }
        }
        
        return result;
    }
    
    template <typename RealType, typename Genotypes>
    RealType probability_haplotype_in_sample(const Haplotype& haplotype, const Genotypes& genotypes,
                                             const SampleGenotypeProbabilities<RealType>& genotype_probabilities)
    {
        RealType result {0};
        
        for (const auto& genotype : genotypes) {
            if (genotype.contains(haplotype)) {
                result += genotype_probabilities.at(genotype);
            }
        }
        
        return result;
    }
    
    template <typename RealType, typename Genotypes>
    RealType probability_genotype_in_samples(const Genotype<Allele>& the_genotype, const Genotypes& genotypes,
                                             const SampleGenotypeProbabilities<RealType>& genotype_probabilities)
    {
        return 0; // TODO
    }
    
    template <typename RealType, typename Haplotypes, typename Genotypes>
    RealType probability_sample_has_genotype(const Genotype<Allele>& the_genotype, const Haplotypes& haplotypes,
                                             const SampleGenotypeProbabilities<RealType>& genotype_probabilities,
                                             const Genotypes& genotypes)
    {
        return 0; // TODO
    }
    
    template <typename RealType, typename Haplotypes>
    RealType probability_allele_in_samples(const Allele& the_allele, const Haplotypes& haplotypes,
                                           const HaplotypePseudoCounts<RealType>& haplotype_pseudo_counts)
    {
        RealType result {0};
        
        for (const auto& haplotype : haplotypes) {
            if (haplotype.contains(the_allele)) {
                result += haplotype_probability(haplotype, haplotype_pseudo_counts);
            }
        }
        
        return result;
    }
    
    template <typename RealType, typename Haplotypes, typename Genotypes>
    RealType probability_allele_in_sample(const Allele& the_allele, const Haplotypes& haplotypes,
                                          const SampleGenotypeProbabilities<RealType>& genotype_probabilities,
                                          const Genotypes& genotypes, RealType zero_epsilon=1e-20)
    {
        std::vector<std::reference_wrapper<const Haplotype>> containing_haplotypes {};
        containing_haplotypes.reserve(haplotypes.size() / 2);
        
        std::copy_if(std::cbegin(haplotypes), std::cend(haplotypes), std::back_inserter(containing_haplotypes),
                     [&the_allele] (const auto& haplotype) { return haplotype.contains(the_allele); });
        
        RealType result {0};
        
        for (const auto& genotype : genotypes) {
            if (genotype_probabilities.at(genotype) >= zero_epsilon &&
                std::any_of(std::cbegin(genotype), std::cend(genotype), [&containing_haplotypes] (const auto& haplotype) {
                return std::binary_search(containing_haplotypes.cbegin(), containing_haplotypes.cend(), haplotype);
            })) {
                result += genotype_probabilities.at(genotype);
            }
        }
        
        return result;
    }
    
    template <typename RealType>
    RealType haplotype_population_probability(const Haplotype& haplotype,
                                              const HaplotypePseudoCounts<RealType>& haplotype_pseudo_counts)
    {
        return haplotype_pseudo_counts.at(haplotype) / sum_values(haplotype_pseudo_counts);
    }
    
    template <typename RealType>
    RealType posterior_predictive_probability(const std::unordered_map<Haplotype, unsigned>& haplotype_occurences,
                                              const HaplotypePseudoCounts<RealType>& haplotype_pseudo_counts)
    {
        std::vector<RealType> z {}, a {};
        z.reserve(haplotype_pseudo_counts.size());
        a.reserve(haplotype_pseudo_counts.size());
        
        for (const auto haplotype_pseudo_count : haplotype_pseudo_counts) {
            z.push_back(haplotype_occurences.at(haplotype_pseudo_count.first));
            a.push_back(haplotype_pseudo_count.second);
        }
        
        return dirichlet_multinomial<RealType>(z, a);
    }
    
    template <typename RealType>
    RealType posterior_predictive_probability(const Genotype<Haplotype>& genotype,
                                              const HaplotypePseudoCounts<RealType>& haplotype_pseudo_counts)
    {
        std::vector<RealType> z {}, a {};
        z.reserve(haplotype_pseudo_counts.size());
        a.reserve(haplotype_pseudo_counts.size());
        
        for (const auto haplotype_pseudo_count : haplotype_pseudo_counts) {
            z.push_back(genotype.num_occurences(haplotype_pseudo_count.first));
            a.push_back(haplotype_pseudo_count.second);
        }
        
        return dirichlet_multinomial<RealType>(z, a);
    }
    
    template <typename RealType>
    RealType median_haplotype_population_probability(const HaplotypePseudoCounts<RealType>& haplotype_pseudo_counts)
    {
        std::vector<RealType> values(haplotype_pseudo_counts.size());
        
        std::transform(std::cbegin(haplotype_pseudo_counts), std::cend(haplotype_pseudo_counts),
                       values.begin(), [] (const auto& pair) {
                           return pair.second;
                       });
        
        std::size_t n {values.size() / 2};
        std::nth_element(values.begin(), values.begin() + n, values.end());
        
        std::cout << values[n] << std::endl;
        
        return (values.size() % 2 == 0) ?
            (values[n] + values[n + 1]) / (2 * std::accumulate(values.cbegin(), values.cend(), 0))
        :
            values[n] / std::accumulate(values.cbegin(), values.cend(), 0);
    }
    
} // end namespace BayesianGenotypeModel

} // end namespace Octopus

#endif
