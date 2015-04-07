//
//  variational_bayes_genotype_model.cpp
//  Octopus
//
//  Created by Daniel Cooke on 01/04/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "variational_bayes_genotype_model.h"

#include <cmath> // std::exp, std::log
#include <boost/math/special_functions/digamma.hpp>

#include "aligned_read.h"
#include "maths.h"
#include "pair_hmm.h"

#include <iostream> // TEST

VariationalBayesGenotypeModel::VariationalBayesGenotypeModel(ReadModel& read_model, unsigned ploidy)
:
ploidy_ {ploidy},
read_model_ {read_model}
{}

double VariationalBayesGenotypeModel::log_expected_genotype_probability(const Genotype& genotype,
                                                                        const HaplotypePseudoCounts& haplotype_pseudo_counts)
{
    // These cases are just for optimisation; they are functionally equivalent
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

double VariationalBayesGenotypeModel::log_rho(const Genotype& genotype,
                                              const HaplotypePseudoCounts& haplotype_pseudo_counts,
                                              const SampleReads& reads, unsigned sample)
{
    return log_expected_genotype_probability(genotype, haplotype_pseudo_counts) +
            read_model_.log_probability(reads, genotype, sample);
}

double VariationalBayesGenotypeModel::genotype_responsability(const Genotype& genotype,
                                                              const SampleReads& reads,
                                                              const HaplotypePseudoCounts& haplotype_pseudo_counts,
                                                              unsigned sample,
                                                              const Genotypes& all_genotypes)
{
    double log_rho_genotype = log_rho(genotype, haplotype_pseudo_counts, reads, sample);
    
    std::vector<double> log_rho_genotypes {};
    
    for (const auto& g : all_genotypes) {
        log_rho_genotypes.push_back(log_rho(g, haplotype_pseudo_counts, reads, sample));
    }
    
    double log_sum_rho_genotypes = log_sum_exp<double>(log_rho_genotypes.cbegin(), log_rho_genotypes.cend());
    
    return std::exp(log_rho_genotype - log_sum_rho_genotypes);
}

double VariationalBayesGenotypeModel::expected_haplotype_count(const Haplotype& haplotype,
                                                               const SampleGenotypeResponsabilities& sample_genotype_responsabilities)
{
    double result {0};
    
    for (const auto& genotype_responsability : sample_genotype_responsabilities) {
        result += genotype_responsability.second * genotype_responsability.first.num_occurences(haplotype);
    }
    
    return result;
}

double VariationalBayesGenotypeModel::posterior_haplotype_pseudo_count(const Haplotype& haplotype,
                                                                       double prior_pseudo_count,
                                                                       const GenotypeResponsabilities& genotype_responsabilities)
{
    double result {prior_pseudo_count};
    
    for (const auto& sample_genotype_responsabilities : genotype_responsabilities) {
        result += expected_haplotype_count(haplotype, sample_genotype_responsabilities);
    }
    
    return result;
}

double VariationalBayesGenotypeModel::posterior_predictive_probability(const std::unordered_map<Haplotype, unsigned>& haplotype_counts,
                                                                       const HaplotypePseudoCounts& haplotype_pseudo_count) const
{
    std::vector<double> z {}, a {};
    
    for (const auto& p : haplotype_counts) { z.push_back(p.second); };
    for (const auto& p : haplotype_pseudo_count) { a.push_back(p.second); };
    
    return dirichlet_multinomial<double>(z, a);
}

double VariationalBayesGenotypeModel::allele_posterior_probability(const Variant& variant,
                                                                   const Haplotypes& haplotypes,
                                                                   const SampleGenotypeResponsabilities& sample_genotype_responsabilities,
                                                                   const Genotypes& genotypes) const
{
    double result {};
    
    for (const auto& haplotype : haplotypes) {
        if (contains(haplotype, variant)) {
            for (const auto& genotype : genotypes) {
                if (genotype.contains(haplotype)) {
                    result += sample_genotype_responsabilities.at(genotype);
                }
            }
        }
    }
    
    return result;
}

double VariationalBayesGenotypeModel::log_expected_genotype_probability_haploid(const Genotype& genotype,
                                                                                const HaplotypePseudoCounts& haplotype_pseudo_counts) const
{
    const static double ln_1 = std::log(1);
    
    return boost::math::digamma<double>(haplotype_pseudo_counts.at(genotype.at(0)))
            - boost::math::digamma<double>(pseudo_count_sum(haplotype_pseudo_counts)) + ln_1;
}

double VariationalBayesGenotypeModel::log_expected_genotype_probability_diploid(const Genotype& genotype,
                                                                                const HaplotypePseudoCounts& haplotype_pseudo_counts) const
{
    const static double ln_1 = std::log(1);
    const static double ln_2 = std::log(2);
    
    if (genotype.is_homozygous()) {
        return 2 * boost::math::digamma<double>(haplotype_pseudo_counts.at(genotype.at(0))) - 2 *
                boost::math::digamma<double>(pseudo_count_sum(haplotype_pseudo_counts)) + ln_1;
    } else {
        return boost::math::digamma<double>(haplotype_pseudo_counts.at(genotype.at(0))) +
                boost::math::digamma<double>(haplotype_pseudo_counts.at(genotype.at(1))) -
                2 * boost::math::digamma<double>(pseudo_count_sum(haplotype_pseudo_counts)) + ln_2;
    }
}

double VariationalBayesGenotypeModel::log_expected_genotype_probability_triploid(const Genotype& genotype,
                                                                                 const HaplotypePseudoCounts& haplotype_pseudo_counts) const
{
    //TODO
    return 0;
}

double VariationalBayesGenotypeModel::log_expected_genotype_probability_polyploid(const Genotype& genotype,
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
        for (unsigned s {}; s < num_samples; ++s) {
            for (const auto& genotype : the_genotypes) {
                responsabilities[s][genotype] = the_model.genotype_responsability(genotype, the_reads.at(s),
                                                                               posterior_pseudo_counts, s,
                                                                               the_genotypes);
            }
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
