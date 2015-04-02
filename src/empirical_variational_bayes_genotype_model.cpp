//
//  empirical_variational_bayes_genotype_model.cpp
//  Octopus
//
//  Created by Daniel Cooke on 01/04/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "empirical_variational_bayes_genotype_model.h"

#include <cmath>
#include <boost/math/special_functions/digamma.hpp>

#include "aligned_read.h"
#include "maths.h"
#include "pair_hmm.h"

#include <iostream> //TEST

EmpiricalVariationalBayesGenotypeModel::EmpiricalVariationalBayesGenotypeModel(ReadModel& read_model,
                                                                               unsigned ploidy)
:
ploidy_ {ploidy},
read_model_ {read_model}
{}

// E_pi [p(genotype | pi)] = sum {haplotype in genotype} (u(haplotype) * y(a_haplotype)) - ploidy * y(a_all) + ln C
// where u(haplotype) is the number of times the haplotype occurs in the genotypes
// a_haplotype is the pseudo count for the haplotype
// a_all is the sum of the pseudo counts for all haplotypes
// y is the digamma function
// and C is the constant multinomial coefficient described in documentation
double EmpiricalVariationalBayesGenotypeModel::log_expected_genotype_probability(const Genotype& genotype,
                                                                                 const HaplotypePseudoCounts& haplotype_pseudo_counts)
{
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

double EmpiricalVariationalBayesGenotypeModel::log_rho(const Genotype& genotype,
                                                       const HaplotypePseudoCounts& haplotype_pseudo_counts,
                                                       const SampleReads& reads, unsigned sample)
{
    return log_expected_genotype_probability(genotype, haplotype_pseudo_counts) +
            read_model_.log_probability(reads, genotype, sample);
}

double EmpiricalVariationalBayesGenotypeModel::genotype_responsability(const Genotype& genotype,
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

double EmpiricalVariationalBayesGenotypeModel::expected_haplotype_count(const Haplotype& haplotype,
                                                                        const SampleGenotypeResponsabilities& sample_genotype_responsabilities)
{
    double result {0};
    
    for (const auto& genotype_responsability : sample_genotype_responsabilities) {
        result += genotype_responsability.second * genotype_responsability.first.num_occurences(haplotype);
    }
    
    return result;
}

double EmpiricalVariationalBayesGenotypeModel::posterior_haplotype_pseudo_count(const Haplotype& haplotype,
                                                                                double prior_pseudo_count,
                                                                                const GenotypeResponsabilities& genotype_responsabilities)
{
    double result {prior_pseudo_count};
    
    for (const auto& sample_genotype_responsabilities : genotype_responsabilities) {
        result += expected_haplotype_count(haplotype, sample_genotype_responsabilities);
    }
    
    return result;
}

double EmpiricalVariationalBayesGenotypeModel::log_expected_genotype_probability_haploid(const Genotype& genotype,
                                                 const HaplotypePseudoCounts& haplotype_pseudo_counts) const
{
    static double ln_1 = std::log(1);
    return boost::math::digamma<double>(haplotype_pseudo_counts.at(genotype.at(0)))
            - boost::math::digamma<double>(pseudo_count_sum(haplotype_pseudo_counts)) + ln_1;
}

double EmpiricalVariationalBayesGenotypeModel::log_expected_genotype_probability_diploid(const Genotype& genotype,
                                                 const HaplotypePseudoCounts& haplotype_pseudo_counts) const
{
    static double ln_1 = std::log(1);
    static double ln_2 = std::log(2);
    
    if (genotype.is_homozygous()) {
        return 2 * boost::math::digamma<double>(haplotype_pseudo_counts.at(genotype.at(0))) - 2 *
                boost::math::digamma<double>(pseudo_count_sum(haplotype_pseudo_counts)) + ln_1;
    } else {
        return boost::math::digamma<double>(haplotype_pseudo_counts.at(genotype.at(0))) +
                boost::math::digamma<double>(haplotype_pseudo_counts.at(genotype.at(1))) -
                ploidy_ * boost::math::digamma<double>(pseudo_count_sum(haplotype_pseudo_counts)) + ln_2;
    }
}

double EmpiricalVariationalBayesGenotypeModel::log_expected_genotype_probability_triploid(const Genotype& genotype,
                                                  const HaplotypePseudoCounts& haplotype_pseudo_counts) const
{
    //TODO
    return 0;
}

double EmpiricalVariationalBayesGenotypeModel::log_expected_genotype_probability_polyploid(const Genotype& genotype,
                                                   const HaplotypePseudoCounts& haplotype_pseudo_counts) const
{
    //TODO
    return 0;
}

unsigned EmpiricalVariationalBayesGenotypeModel::pseudo_count_sum(const HaplotypePseudoCounts& haplotype_pseudo_counts) const
{
    unsigned result {0};
    for (const auto& h : haplotype_pseudo_counts) {
        result += h.second;
    }
    return result;
}
