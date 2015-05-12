//
//  variational_bayes_genotype_model.h.h
//  Octopus
//
//  Created by Daniel Cooke on 01/04/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__variational_bayes_genotype_model__
#define __Octopus__variational_bayes_genotype_model__

#include <vector>
#include <string>
#include <unordered_map>
#include <algorithm>  // std::transform
#include <cmath>      // std::exp, std::log

#include "common.h"
#include "haplotype.h"
#include "genotype.h"
#include "read_model.h"
#include "maths.h"
#include "aligned_read.h"
#include "pair_hmm.h"
#include "bayesian_genotype_model.h"

namespace BayesianGenotypeModel
{
    template <typename ForwardIterator>
    using SampleReadRange = std::pair<ForwardIterator, ForwardIterator>;
    
    template <typename SampleIdType, typename ForwardIterator>
    using ReadRanges = std::unordered_map<SampleIdType, SampleReadRange<ForwardIterator>>;
    
    class VariationalBayesGenotypeModel
    {
    public:
        using RealType     = Octopus::ProbabilityType;
        using SampleIdType = Octopus::SampleIdType;
        
        VariationalBayesGenotypeModel() = delete;
        explicit VariationalBayesGenotypeModel(ReadModel& read_model, unsigned ploidy);
        ~VariationalBayesGenotypeModel() = default;
        
        VariationalBayesGenotypeModel(const VariationalBayesGenotypeModel&)            = default;
        VariationalBayesGenotypeModel& operator=(const VariationalBayesGenotypeModel&) = default;
        VariationalBayesGenotypeModel(VariationalBayesGenotypeModel&&)                 = default;
        VariationalBayesGenotypeModel& operator=(VariationalBayesGenotypeModel&&)      = default;
        
        RealType log_expected_genotype_probability(const Genotype& genotype,
                                                   const HaplotypePseudoCounts<RealType>& haplotype_pseudo_counts);
        
        template <typename ForwardIterator>
        RealType log_rho(const Genotype& genotype, const HaplotypePseudoCounts<RealType>& haplotype_pseudo_counts,
                         ForwardIterator first_read, ForwardIterator last_read, SampleIdType sample);
        
        template <typename ForwardIterator, typename Container>
        RealType genotype_posterior(const Genotype& genotype, ForwardIterator first_read, ForwardIterator last_read,
                                    const HaplotypePseudoCounts<RealType>& haplotype_pseudo_counts,
                                    const Container& genotypes, SampleIdType sample);
        
        template <typename ForwardIterator, typename Container>
        SampleGenotypeProbabilities<RealType> genotype_posteriors(const Container& genotypes,
                                                               ForwardIterator first_read, ForwardIterator last_read,
                                                               const HaplotypePseudoCounts<RealType>& haplotype_pseudo_counts,
                                                               SampleIdType sample);
        
        RealType expected_haplotype_count(const Haplotype& haplotype,
                                          const SampleGenotypeProbabilities<RealType>& sample_genotype_responsabilities) const;
        
        RealType posterior_haplotype_pseudo_count(const Haplotype& haplotype, RealType prior_pseudo_count,
                                                  const GenotypeProbabilities<SampleIdType, RealType>& genotype_responsabilities) const;
        
        // This is just a slight optimisation of the other posterior_haplotype_pseudo_count
        template <typename Container>
        RealType posterior_haplotype_pseudo_count(const Haplotype& haplotype, RealType prior_pseudo_count,
                                                  const GenotypeProbabilities<SampleIdType, RealType>& genotype_responsabilities,
                                                  const Container& genotypes) const;
        
        void clear_cache();
        
    private:
        unsigned ploidy_;
        ReadModel& read_model_;
        
        // These are just for optimisation
        RealType log_expected_genotype_probability_haploid(const Genotype& genotype,
                                                           const HaplotypePseudoCounts<RealType>& haplotype_pseudo_counts) const;
        RealType log_expected_genotype_probability_diploid(const Genotype& genotype,
                                                           const HaplotypePseudoCounts<RealType>& haplotype_pseudo_counts) const;
        RealType log_expected_genotype_probability_triploid(const Genotype& genotype,
                                                            const HaplotypePseudoCounts<RealType>& haplotype_pseudo_counts) const;
        RealType log_expected_genotype_probability_polyploid(const Genotype& genotype,
                                                             const HaplotypePseudoCounts<RealType>& haplotype_pseudo_counts) const;
    };
    
    template <typename ForwardIterator>
    inline
    VariationalBayesGenotypeModel::RealType
    VariationalBayesGenotypeModel::log_rho(const Genotype& genotype,
                                           const HaplotypePseudoCounts<RealType>& haplotype_pseudo_counts,
                                           ForwardIterator first_read, ForwardIterator last_read,
                                           SampleIdType sample)
    {
        return log_expected_genotype_probability(genotype, haplotype_pseudo_counts) +
        read_model_.log_probability(first_read, last_read, genotype, sample);
    }
    
    template <typename ForwardIterator, typename Container>
    VariationalBayesGenotypeModel::RealType
    VariationalBayesGenotypeModel::genotype_posterior(const Genotype& genotype,
                                                      ForwardIterator first_read, ForwardIterator last_read,
                                                      const HaplotypePseudoCounts<RealType>& haplotype_pseudo_counts,
                                                      const Container& genotypes, SampleIdType sample)
    {
        RealType log_rho_genotype = log_rho(genotype, haplotype_pseudo_counts, first_read, last_read, sample);
        
        std::vector<RealType> log_rho_genotypes (genotypes.size());
        
        std::transform(std::cbegin(genotypes), std::cend(genotypes), log_rho_genotypes.begin(),
                       [this, &haplotype_pseudo_counts, first_read, last_read, sample] (const auto& genotype) {
                           return log_rho(genotype, haplotype_pseudo_counts, first_read, last_read, sample);
                       });
        
        RealType log_sum_rho_genotypes = log_sum_exp<RealType>(log_rho_genotypes.cbegin(), log_rho_genotypes.cend());
        
        return std::exp(log_rho_genotype - log_sum_rho_genotypes);
    }
    
    template <typename ForwardIterator, typename Container>
    SampleGenotypeProbabilities<VariationalBayesGenotypeModel::RealType>
    VariationalBayesGenotypeModel::genotype_posteriors(const Container& genotypes,
                                                       ForwardIterator first_read, ForwardIterator last_read,
                                                       const HaplotypePseudoCounts<RealType>& haplotype_pseudo_counts,
                                                       SampleIdType sample)
    {
        std::vector<RealType> log_rho_genotypes (genotypes.size());
        
        std::transform(std::cbegin(genotypes), std::cend(genotypes), log_rho_genotypes.begin(),
                       [this, &haplotype_pseudo_counts, first_read, last_read, sample] (const auto& genotype) {
                           return log_rho(genotype, haplotype_pseudo_counts, first_read, last_read, sample);
                       });
        
        RealType log_sum_rho_genotypes = log_sum_exp<RealType>(log_rho_genotypes.cbegin(), log_rho_genotypes.cend());
        
        SampleGenotypeProbabilities<RealType> result {};
        result.reserve(genotypes.size());
        
        for (unsigned i {}; i < genotypes.size(); ++i) {
            result[genotypes.at(i)] = std::exp(log_rho_genotypes.at(i) - log_sum_rho_genotypes);
        }
        
        return result;
    }
    
    template <typename Container>
    VariationalBayesGenotypeModel::RealType
    VariationalBayesGenotypeModel::posterior_haplotype_pseudo_count(const Haplotype& haplotype,
                                                                    RealType prior_pseudo_count,
                                                                    const GenotypeProbabilities<SampleIdType, RealType>& genotype_responsabilities,
                                                                    const Container& genotypes) const
    {
        RealType result {prior_pseudo_count};
        
        unsigned num_occurences {};
        RealType responsability_sum {};
        
        for (const auto& genotype : genotypes) {
            num_occurences = genotype.num_occurences(haplotype);
            
            if (num_occurences > 0) {
                for (const auto& sample_genotype_responsabilities : genotype_responsabilities) {
                    responsability_sum += sample_genotype_responsabilities.second.at(genotype);
                }
                
                result += num_occurences * responsability_sum;
                responsability_sum = 0;
            }
        }
        
        return result;
    }
    
    // Non-member functions
    
    template <typename RealType, typename SampleIdType, typename ForwardIterator, typename Container>
    Latents<SampleIdType, RealType>
    update_latents(VariationalBayesGenotypeModel& the_model, const Container& the_genotypes,
                   const HaplotypePseudoCounts<RealType>& haplotype_prior_counts,
                   const ReadRanges<SampleIdType, ForwardIterator>& the_reads, unsigned max_num_iterations)
    {
        GenotypeProbabilities<SampleIdType, RealType> genotype_posteriors(the_reads.size());
        HaplotypePseudoCounts<RealType> haplotype_posterior_counts {haplotype_prior_counts};
        
        for (unsigned i {}; i < max_num_iterations; ++i) {
            for (const auto& sample : the_reads) {
                genotype_posteriors[sample.first] =
                the_model.genotype_posteriors(the_genotypes, sample.second.first,
                                              sample.second.second, haplotype_posterior_counts,
                                              sample.first);
            }
            
            for (const auto& haplotype_prior : haplotype_prior_counts) {
                haplotype_posterior_counts[haplotype_prior.first] =
                the_model.posterior_haplotype_pseudo_count(haplotype_prior.first, haplotype_prior.second,
                                                           genotype_posteriors, the_genotypes);
            }
        }
        
        return Latents<SampleIdType, RealType> {std::move(haplotype_posterior_counts),
                                                std::move(genotype_posteriors)};
    }
    
} // end namespace BayesianGenotypeModel

#endif /* defined(__Octopus__variational_bayes_genotype_model.h__) */
