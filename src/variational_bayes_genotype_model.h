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

#include "haplotype.h"
#include "genotype.h"
#include "read_model.h"
#include "maths.h"

class AlignedRead;

class VariationalBayesGenotypeModel
{
public:
    using RealType                       = double;
    using Haplotypes                     = std::vector<Haplotype>;
    using Genotypes                      = std::vector<Genotype>;
    using HaplotypePseudoCounts          = std::unordered_map<Haplotype, RealType>;
    using SampleGenotypeResponsabilities = std::unordered_map<Genotype, RealType>;
    using SampleIdType                   = std::string;
    using GenotypeResponsabilities       = std::unordered_map<SampleIdType, SampleGenotypeResponsabilities>;
    
    VariationalBayesGenotypeModel() = delete;
    explicit VariationalBayesGenotypeModel(ReadModel& read_model, unsigned ploidy,
                                           RealType zero_epsilon=1e-20);
    ~VariationalBayesGenotypeModel() = default;
    
    VariationalBayesGenotypeModel(const VariationalBayesGenotypeModel&)            = default;
    VariationalBayesGenotypeModel& operator=(const VariationalBayesGenotypeModel&) = default;
    VariationalBayesGenotypeModel(VariationalBayesGenotypeModel&&)                 = default;
    VariationalBayesGenotypeModel& operator=(VariationalBayesGenotypeModel&&)      = default;
    
    RealType log_expected_genotype_probability(const Genotype& genotype,
                                               const HaplotypePseudoCounts& haplotype_pseudo_counts);
    
    template <typename ForwardIterator>
    RealType log_rho(const Genotype& genotype, const HaplotypePseudoCounts& haplotype_pseudo_counts,
                     ForwardIterator first_read, ForwardIterator last_read, SampleIdType sample);
    
    template <typename ForwardIterator>
    RealType genotype_responsability(const Genotype& genotype, ForwardIterator first_read, ForwardIterator last_read,
                                     const HaplotypePseudoCounts& haplotype_pseudo_counts,
                                     const Genotypes& genotypes, SampleIdType sample);
    
    template <typename ForwardIterator>
    SampleGenotypeResponsabilities genotype_responsabilities(const Genotypes& genotypes,
                                                             ForwardIterator first_read, ForwardIterator last_read,
                                                             const HaplotypePseudoCounts& haplotype_pseudo_counts,
                                                             SampleIdType sample);
    
    RealType expected_haplotype_count(const Haplotype& haplotype,
                                      const SampleGenotypeResponsabilities& sample_genotype_responsabilities) const;
    
    RealType posterior_haplotype_pseudo_count(const Haplotype& haplotype, RealType prior_pseudo_count,
                                              const GenotypeResponsabilities& genotype_responsabilities) const;
    
    // This is just a slight optimisation of the other posterior_haplotype_pseudo_count
    RealType posterior_haplotype_pseudo_count(const Haplotype& haplotype, RealType prior_pseudo_count,
                                              const GenotypeResponsabilities& genotype_responsabilities,
                                              const Genotypes& genotypes) const;
    
    RealType posterior_haplotype_probability(const Haplotype& haplotype,
                                             const HaplotypePseudoCounts& posterior_haplotype_pseudo_counts) const;
    
    RealType posterior_predictive_probability(const std::unordered_map<Haplotype, unsigned>& haplotype_counts,
                                              const HaplotypePseudoCounts& haplotype_pseudo_counts) const;
    
    RealType posterior_predictive_probability(const Genotype& genotype,
                                              const HaplotypePseudoCounts& haplotype_pseudo_counts) const;
    
    RealType posterior_probability_haplotype_in_samples(const Haplotype& haplotype,
                                                        const Genotypes& genotypes,
                                                        const HaplotypePseudoCounts& posterior_haplotype_pseudo_counts) const;
    
    RealType posterior_probability_haplotype_in_sample(const Haplotype& haplotype,
                                                       const Genotypes& genotypes,
                                                       const SampleGenotypeResponsabilities& genotype_responsabilities) const;
    
    RealType posterior_probability_allele_in_samples(const Allele& the_allele,
                                                     const Haplotypes& haplotypes,
                                                     const HaplotypePseudoCounts& posterior_haplotype_pseudo_counts) const;
    
    RealType posterior_probability_allele_in_sample(const Allele& the_allele,
                                                    const Haplotypes& haplotypes,
                                                    const SampleGenotypeResponsabilities& sample_genotype_responsabilities,
                                                    const Genotypes& genotypes) const;
    void clear_cache();
    
private:
    unsigned ploidy_;
    ReadModel& read_model_;
    RealType zero_epsilon_;
    
    // These are just for optimisation
    RealType log_expected_genotype_probability_haploid(const Genotype& genotype,
                                                       const HaplotypePseudoCounts& haplotype_pseudo_counts) const;
    RealType log_expected_genotype_probability_diploid(const Genotype& genotype,
                                                       const HaplotypePseudoCounts& haplotype_pseudo_counts) const;
    RealType log_expected_genotype_probability_triploid(const Genotype& genotype,
                                                        const HaplotypePseudoCounts& haplotype_pseudo_counts) const;
    RealType log_expected_genotype_probability_polyploid(const Genotype& genotype,
                                                         const HaplotypePseudoCounts& haplotype_pseudo_counts) const;
};

template <typename ForwardIterator>
inline
VariationalBayesGenotypeModel::RealType
VariationalBayesGenotypeModel::log_rho(const Genotype& genotype,
                                       const HaplotypePseudoCounts& haplotype_pseudo_counts,
                                       ForwardIterator first_read, ForwardIterator last_read,
                                       SampleIdType sample)
{
    return log_expected_genotype_probability(genotype, haplotype_pseudo_counts) +
    read_model_.log_probability(first_read, last_read, genotype, sample);
}

template <typename ForwardIterator>
VariationalBayesGenotypeModel::RealType
VariationalBayesGenotypeModel::genotype_responsability(const Genotype& genotype,
                                                       ForwardIterator first_read, ForwardIterator last_read,
                                                       const HaplotypePseudoCounts& haplotype_pseudo_counts,
                                                       const Genotypes& genotypes, SampleIdType sample)
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

template <typename ForwardIterator>
VariationalBayesGenotypeModel::SampleGenotypeResponsabilities
VariationalBayesGenotypeModel::genotype_responsabilities(const Genotypes& genotypes,
                                                         ForwardIterator first_read, ForwardIterator last_read,
                                                         const HaplotypePseudoCounts& haplotype_pseudo_counts,
                                                         SampleIdType sample)
{
    std::vector<RealType> log_rho_genotypes (genotypes.size());
    
    std::transform(std::cbegin(genotypes), std::cend(genotypes), log_rho_genotypes.begin(),
                   [this, &haplotype_pseudo_counts, first_read, last_read, sample] (const auto& genotype) {
                       return log_rho(genotype, haplotype_pseudo_counts, first_read, last_read, sample);
                   });
    
    RealType log_sum_rho_genotypes = log_sum_exp<RealType>(log_rho_genotypes.cbegin(), log_rho_genotypes.cend());
    
    SampleGenotypeResponsabilities result {};
    result.reserve(genotypes.size());
    
    for (unsigned i {}; i < genotypes.size(); ++i) {
        result[genotypes.at(i)] = std::exp(log_rho_genotypes.at(i) - log_sum_rho_genotypes);
    }
    
    return result;
}

template <typename T, typename RealType>
RealType sum(const std::unordered_map<T, RealType>& map) noexcept
{
    RealType result {};
    
    for (const auto& p : map) {
        result += p.second;
    }
    
    return result;
}

using HaplotypePriors = std::unordered_map<Haplotype, VariationalBayesGenotypeModel::RealType>;

VariationalBayesGenotypeModel::HaplotypePseudoCounts
get_haplotype_prior_pseudo_counts(const HaplotypePriors& the_haplotype_priors,
                                  const Haplotype& the_reference_haplotype,
                                  VariationalBayesGenotypeModel::RealType the_reference_haplotype_pseudo_count);

template <typename ForwardIterator>
using SamplesReads = std::unordered_map<VariationalBayesGenotypeModel::SampleIdType,
                                        std::pair<ForwardIterator, ForwardIterator>>;

using VariationalBayesGenotypeModelLatents = std::pair<VariationalBayesGenotypeModel::GenotypeResponsabilities,
                                                    VariationalBayesGenotypeModel::HaplotypePseudoCounts>;

template <typename ForwardIterator>
VariationalBayesGenotypeModelLatents
update_parameters(VariationalBayesGenotypeModel& the_model,
                  const VariationalBayesGenotypeModel::Genotypes& the_genotypes,
                  const VariationalBayesGenotypeModel::HaplotypePseudoCounts& prior_haplotype_pseudocounts,
                  const SamplesReads<ForwardIterator>& the_reads,
                  unsigned max_num_iterations)
{
    VariationalBayesGenotypeModel::GenotypeResponsabilities responsabilities(the_reads.size());
    VariationalBayesGenotypeModel::HaplotypePseudoCounts posterior_pseudo_counts {prior_haplotype_pseudocounts};
    
    for (unsigned i {}; i < max_num_iterations; ++i) {
        for (const auto& sample : the_reads) {
            responsabilities[sample.first] = the_model.genotype_responsabilities(the_genotypes,
                                                                                 sample.second.first,
                                                                                 sample.second.second,
                                                                                 posterior_pseudo_counts,
                                                                                 sample.first);
        }
        
        for (const auto& haplotype_prior : prior_haplotype_pseudocounts) {
            posterior_pseudo_counts[haplotype_prior.first] = the_model.posterior_haplotype_pseudo_count(haplotype_prior.first,
                                                                                                        haplotype_prior.second,
                                                                                                        responsabilities,
                                                                                                        the_genotypes);
        }
    }
    
    return {responsabilities, posterior_pseudo_counts};
}

#endif /* defined(__Octopus__variational_bayes_genotype_model.h__) */
