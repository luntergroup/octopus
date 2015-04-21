//
//  variational_bayes_genotype_model.cpp
//  Octopus
//
//  Created by Daniel Cooke on 01/04/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "variational_bayes_genotype_model.h"

#include <cmath>      // std::exp, std::log
#include <algorithm>  // std::transform, std::copy_if, std::any_of
#include <functional> // std::reference_wrapper
#include <boost/math/special_functions/digamma.hpp>

#include "aligned_read.h"
#include "maths.h"
#include "pair_hmm.h"

VariationalBayesGenotypeModel::VariationalBayesGenotypeModel(ReadModel& read_model, unsigned ploidy,
                                                             RealType zero_epsilon)
:
ploidy_ {ploidy},
read_model_ {read_model},
zero_epsilon_ {zero_epsilon}
{}

VariationalBayesGenotypeModel::RealType
VariationalBayesGenotypeModel::log_expected_genotype_probability(const Genotype& genotype,
                                                                 const HaplotypePseudoCounts& haplotype_pseudo_counts)
{
    // These cases are just for optimisation
    switch (ploidy_) {
        case 1:
            return log_expected_genotype_probability_haploid(genotype, haplotype_pseudo_counts);
        case 2:
            return log_expected_genotype_probability_diploid(genotype, haplotype_pseudo_counts);
        case 3:
            return log_expected_genotype_probability_triploid(genotype, haplotype_pseudo_counts);
        default:
            return log_expected_genotype_probability_polyploid(genotype, haplotype_pseudo_counts);
    }
}

VariationalBayesGenotypeModel::RealType
VariationalBayesGenotypeModel::log_rho(const Genotype& genotype,
                                       const HaplotypePseudoCounts& haplotype_pseudo_counts,
                                       ReadIterator first, ReadIterator last, unsigned sample)
{
    return log_expected_genotype_probability(genotype, haplotype_pseudo_counts) +
            read_model_.log_probability(first, last, genotype, sample);
}

VariationalBayesGenotypeModel::RealType
VariationalBayesGenotypeModel::genotype_responsability(const Genotype& genotype,
                                                       ReadIterator first, ReadIterator last,
                                                       const HaplotypePseudoCounts& haplotype_pseudo_counts,
                                                       unsigned sample,
                                                       const Genotypes& genotypes)
{
    RealType log_rho_genotype = log_rho(genotype, haplotype_pseudo_counts, first, last, sample);
    
    std::vector<RealType> log_rho_genotypes (genotypes.size());
    
    std::transform(std::cbegin(genotypes), std::cend(genotypes), log_rho_genotypes.begin(),
                   [this, &haplotype_pseudo_counts, first, last, sample] (const auto& genotype) {
                       return log_rho(genotype, haplotype_pseudo_counts, first, last, sample);
                   });
    
    RealType log_sum_rho_genotypes = log_sum_exp<RealType>(log_rho_genotypes.cbegin(), log_rho_genotypes.cend());
    
    return std::exp(log_rho_genotype - log_sum_rho_genotypes);
}

VariationalBayesGenotypeModel::SampleGenotypeResponsabilities
VariationalBayesGenotypeModel::genotype_responsabilities(const Genotypes& genotypes, ReadIterator first,
                                                         ReadIterator last,
                                                         const HaplotypePseudoCounts& haplotype_pseudo_counts,
                                                         unsigned sample)
{
    std::vector<RealType> log_rho_genotypes (genotypes.size());
    
    std::transform(std::cbegin(genotypes), std::cend(genotypes), log_rho_genotypes.begin(),
                   [this, &haplotype_pseudo_counts, first, last, sample] (const auto& genotype) {
                       return log_rho(genotype, haplotype_pseudo_counts, first, last, sample);
                   });
    
    RealType log_sum_rho_genotypes = log_sum_exp<RealType>(log_rho_genotypes.cbegin(), log_rho_genotypes.cend());
    
    SampleGenotypeResponsabilities result {};
    result.reserve(genotypes.size());
    
    for (unsigned i {}; i < genotypes.size(); ++i) {
        result[genotypes.at(i)] = std::exp(log_rho_genotypes.at(i) - log_sum_rho_genotypes);
    }
    
    return result;
}

VariationalBayesGenotypeModel::RealType
VariationalBayesGenotypeModel::expected_haplotype_count(const Haplotype& haplotype,
                                                        const SampleGenotypeResponsabilities& sample_genotype_responsabilities)
{
    RealType result {0};
    
    for (const auto& genotype_responsability : sample_genotype_responsabilities) {
        result += genotype_responsability.second * genotype_responsability.first.num_occurences(haplotype);
    }
    
    return result;
}

VariationalBayesGenotypeModel::RealType
VariationalBayesGenotypeModel::posterior_haplotype_pseudo_count(const Haplotype& haplotype,
                                                                RealType prior_pseudo_count,
                                                                const GenotypeResponsabilities& genotype_responsabilities)
{
    RealType result {prior_pseudo_count};
    
    for (const auto& sample_genotype_responsabilities : genotype_responsabilities) {
        result += expected_haplotype_count(haplotype, sample_genotype_responsabilities);
    }
    
    return result;
}

VariationalBayesGenotypeModel::RealType
VariationalBayesGenotypeModel::posterior_haplotype_probability(const Haplotype& haplotype,
                                                               const HaplotypePseudoCounts& posterior_haplotype_pseudo_counts) const
{
    return posterior_haplotype_pseudo_counts.at(haplotype) / pseudo_count_sum(posterior_haplotype_pseudo_counts);
}

VariationalBayesGenotypeModel::RealType
VariationalBayesGenotypeModel::posterior_predictive_probability(const std::unordered_map<Haplotype, unsigned>& haplotype_counts,
                                                                const HaplotypePseudoCounts& haplotype_pseudo_counts) const
{
    std::vector<RealType> z {}, a {};
    z.reserve(haplotype_pseudo_counts.size());
    a.reserve(haplotype_pseudo_counts.size());
    
    for (const auto haplotype_pseudo_count : haplotype_pseudo_counts) {
        z.push_back(haplotype_counts.at(haplotype_pseudo_count.first));
        a.push_back(haplotype_pseudo_count.second);
    }
    
    return dirichlet_multinomial<RealType>(z, a);
}

VariationalBayesGenotypeModel::RealType
VariationalBayesGenotypeModel::posterior_predictive_probability(const Genotype& genotype,
                                                                const HaplotypePseudoCounts& haplotype_pseudo_counts) const
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

VariationalBayesGenotypeModel::RealType
VariationalBayesGenotypeModel::posterior_probability_haplotype_in_samples(const Haplotype& haplotype,
                                                                          const Genotypes& all_genotypes,
                                                                          const HaplotypePseudoCounts& posterior_haplotype_pseudo_counts) const
{
    RealType result {0};
    
    for (const auto& genotype : all_genotypes) {
        if (genotype.contains(haplotype)) {
            result += posterior_predictive_probability(genotype, posterior_haplotype_pseudo_counts);
        }
    }
    
    return result;
}

VariationalBayesGenotypeModel::RealType
VariationalBayesGenotypeModel::posterior_probability_haplotype_in_sample(const Haplotype& haplotype,
                                                                         const Genotypes& all_genotypes,
                                                                         const SampleGenotypeResponsabilities& genotype_responsabilities) const
{
    RealType result {0};
    
    for (const auto& genotype : all_genotypes) {
        if (genotype.contains(haplotype)) {
            result += genotype_responsabilities.at(genotype);
        }
    }
    
    return result;
}

VariationalBayesGenotypeModel::RealType
VariationalBayesGenotypeModel::posterior_probability_allele_in_samples(const Allele& the_allele,
                                                                       const Haplotypes& haplotypes,
                                                                       const HaplotypePseudoCounts& posterior_haplotype_pseudo_counts) const
{
    RealType result {0};
    
    for (const auto& haplotype : haplotypes) {
        if (haplotype.contains(the_allele)) {
            result += posterior_haplotype_probability(haplotype, posterior_haplotype_pseudo_counts);
        }
    }
    
    return result;
}

VariationalBayesGenotypeModel::RealType
VariationalBayesGenotypeModel::posterior_probability_allele_in_sample(const Allele& the_allele,
                                                                      const Haplotypes& haplotypes,
                                                                      const SampleGenotypeResponsabilities& sample_genotype_responsabilities,
                                                                      const Genotypes& genotypes) const
{
    std::vector<std::reference_wrapper<const Haplotype>> containing_haplotypes {};
    containing_haplotypes.reserve(haplotypes.size() / 2);
    
    std::copy_if(cbegin(haplotypes), std::cend(haplotypes), std::back_inserter(containing_haplotypes),
                 [&the_allele] (const auto& haplotype) { return haplotype.contains(the_allele); });
    
    RealType result {0};
    
    for (const auto& genotype : genotypes) {
        if (sample_genotype_responsabilities.at(genotype) >= zero_epsilon_ &&
            std::any_of(containing_haplotypes.cbegin(), containing_haplotypes.cend(),
                        [&genotype] (const auto& haplotype) { return genotype.contains(haplotype); })) {
            result += sample_genotype_responsabilities.at(genotype);
        }
    }
    
    return result;
}

// Private methods

VariationalBayesGenotypeModel::RealType
VariationalBayesGenotypeModel::log_expected_genotype_probability_haploid(const Genotype& genotype,
                                                                         const HaplotypePseudoCounts& haplotype_pseudo_counts) const
{
    const static RealType ln_1 = std::log(1);
    
    return boost::math::digamma<RealType>(haplotype_pseudo_counts.at(genotype.at(0)))
            - boost::math::digamma<RealType>(pseudo_count_sum(haplotype_pseudo_counts)) + ln_1;
}

VariationalBayesGenotypeModel::RealType
VariationalBayesGenotypeModel::log_expected_genotype_probability_diploid(const Genotype& genotype,
                                                                         const HaplotypePseudoCounts& haplotype_pseudo_counts) const
{
    const static RealType ln_1 = std::log(1);
    const static RealType ln_2 = std::log(2);
    
    return (genotype.is_homozygous()) ?
        2 * boost::math::digamma<RealType>(haplotype_pseudo_counts.at(genotype.at(0))) -
        2 * boost::math::digamma<RealType>(pseudo_count_sum(haplotype_pseudo_counts)) + ln_1
        :
        boost::math::digamma<RealType>(haplotype_pseudo_counts.at(genotype.at(0))) +
        boost::math::digamma<RealType>(haplotype_pseudo_counts.at(genotype.at(1))) -
        2 * boost::math::digamma<RealType>(pseudo_count_sum(haplotype_pseudo_counts)) + ln_2;
}

VariationalBayesGenotypeModel::RealType
VariationalBayesGenotypeModel::log_expected_genotype_probability_triploid(const Genotype& genotype,
                                                                          const HaplotypePseudoCounts& haplotype_pseudo_counts) const
{
    //TODO
    return 0;
}

VariationalBayesGenotypeModel::RealType
VariationalBayesGenotypeModel::log_expected_genotype_probability_polyploid(const Genotype& genotype,
                                                                           const HaplotypePseudoCounts& haplotype_pseudo_counts) const
{
    //TODO
    return 0;
}

unsigned VariationalBayesGenotypeModel::pseudo_count_sum(const HaplotypePseudoCounts& haplotype_pseudo_counts) const
{
    unsigned result {0};
    
    for (const auto& h : haplotype_pseudo_counts) {
        result += h.second;
    }
    
    return result;
}

GenotypePosteriors update_parameters(VariationalBayesGenotypeModel& the_model,
                                     const VariationalBayesGenotypeModel::Genotypes& the_genotypes,
                                     const VariationalBayesGenotypeModel::HaplotypePseudoCounts& prior_haplotype_pseudocounts,
                                     const SamplesReads& the_reads, unsigned max_num_iterations)
{
    unsigned num_samples {static_cast<unsigned>(the_reads.size())};
    
    VariationalBayesGenotypeModel::GenotypeResponsabilities responsabilities(num_samples);
    VariationalBayesGenotypeModel::HaplotypePseudoCounts posterior_pseudo_counts {prior_haplotype_pseudocounts};
    
    for (unsigned i {}; i < max_num_iterations; ++i) {
        for (unsigned sample {}; sample < num_samples; ++sample) {
            responsabilities[sample] = the_model.genotype_responsabilities(the_genotypes,
                                                                      the_reads.at(sample).first,
                                                                      the_reads.at(sample).second,
                                                                      posterior_pseudo_counts, sample);
        }
        
        for (const auto& haplotype_prior_pair : prior_haplotype_pseudocounts) {
            const auto& haplotype = haplotype_prior_pair.first;
            posterior_pseudo_counts[haplotype] = the_model.posterior_haplotype_pseudo_count(haplotype,
                                                                                            haplotype_prior_pair.second,
                                                                                            responsabilities);
        }
    }
    
    return {responsabilities, posterior_pseudo_counts};
}
