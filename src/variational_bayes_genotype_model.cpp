//
//  variational_bayes_genotype_model.cpp
//  Octopus
//
//  Created by Daniel Cooke on 01/04/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "variational_bayes_genotype_model.h"

#include <cmath>      // std::exp, std::log
#include <algorithm>  // std::transform, std::copy_if, std::any_of, std::binary_search
#include <functional> // std::reference_wrapper
#include <boost/math/special_functions/digamma.hpp>

#include "aligned_read.h"
#include "pair_hmm.h"

#include <iostream> // TEST

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
VariationalBayesGenotypeModel::expected_haplotype_count(const Haplotype& haplotype,
                                                        const SampleGenotypeResponsabilities& sample_genotype_responsabilities) const
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
                                                                const GenotypeResponsabilities& genotype_responsabilities) const
{
    RealType result {prior_pseudo_count};
    
    for (const auto& sample_genotype_responsabilities : genotype_responsabilities) {
        result += expected_haplotype_count(haplotype, sample_genotype_responsabilities.second);
    }
    
    return result;
}

VariationalBayesGenotypeModel::RealType
VariationalBayesGenotypeModel::posterior_haplotype_pseudo_count(const Haplotype& haplotype,
                                                                RealType prior_pseudo_count,
                                                                const GenotypeResponsabilities& genotype_responsabilities,
                                                                const Genotypes& genotypes) const
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

VariationalBayesGenotypeModel::RealType
VariationalBayesGenotypeModel::posterior_haplotype_probability(const Haplotype& haplotype,
                                                               const HaplotypePseudoCounts& posterior_haplotype_pseudo_counts) const
{
    return posterior_haplotype_pseudo_counts.at(haplotype) / sum(posterior_haplotype_pseudo_counts);
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

// It is required that haplotypes is sorted
VariationalBayesGenotypeModel::RealType
VariationalBayesGenotypeModel::posterior_probability_allele_in_sample(const Allele& the_allele,
                                                                      const Haplotypes& haplotypes,
                                                                      const SampleGenotypeResponsabilities& sample_genotype_responsabilities,
                                                                      const Genotypes& genotypes) const
{
    std::vector<std::reference_wrapper<const Haplotype>> containing_haplotypes {};
    containing_haplotypes.reserve(haplotypes.size() / 2);
    
    std::copy_if(std::cbegin(haplotypes), std::cend(haplotypes), std::back_inserter(containing_haplotypes),
                 [&the_allele] (const auto& haplotype) { return haplotype.contains(the_allele); });
    
    RealType result {0};
    
    for (const auto& genotype : genotypes) {
        if (sample_genotype_responsabilities.at(genotype) >= zero_epsilon_ &&
            std::any_of(std::cbegin(genotype), std::cend(genotype), [&containing_haplotypes] (const auto& haplotype) {
                return std::binary_search(containing_haplotypes.cbegin(), containing_haplotypes.cend(), haplotype);
            })) {
                result += sample_genotype_responsabilities.at(genotype);
        }
    }
    
    return result;
}

void VariationalBayesGenotypeModel::clear_cache()
{
    read_model_.clear_cache();
}

// Private methods

VariationalBayesGenotypeModel::RealType
VariationalBayesGenotypeModel::log_expected_genotype_probability_haploid(const Genotype& genotype,
                                                                         const HaplotypePseudoCounts& haplotype_pseudo_counts) const
{
    return boost::math::digamma<RealType>(haplotype_pseudo_counts.at(genotype.at(0)))
            - boost::math::digamma<RealType>(sum(haplotype_pseudo_counts));
}

VariationalBayesGenotypeModel::RealType
VariationalBayesGenotypeModel::log_expected_genotype_probability_diploid(const Genotype& genotype,
                                                                         const HaplotypePseudoCounts& haplotype_pseudo_counts) const
{
    const static RealType ln_2 = std::log(2);
    
    return ((genotype.is_homozygous()) ?
        2 * boost::math::digamma<RealType>(haplotype_pseudo_counts.at(genotype.at(0))) :
        ln_2 + boost::math::digamma<RealType>(haplotype_pseudo_counts.at(genotype.at(0))) +
            boost::math::digamma<RealType>(haplotype_pseudo_counts.at(genotype.at(1))))
        - 2 * boost::math::digamma<RealType>(sum(haplotype_pseudo_counts));
}

VariationalBayesGenotypeModel::RealType
VariationalBayesGenotypeModel::log_expected_genotype_probability_triploid(const Genotype& genotype,
                                                                          const HaplotypePseudoCounts& haplotype_pseudo_counts) const
{
    const static RealType ln_3 = std::log(3);
    const static RealType ln_6 = std::log(6);
    
    auto k = 3 * boost::math::digamma<RealType>(sum(haplotype_pseudo_counts));
    
    if (genotype.is_homozygous()) {
        return 3 * boost::math::digamma<RealType>(haplotype_pseudo_counts.at(genotype.at(0))) - k;
    } else if (genotype.num_occurences(genotype.at(0)) == 2) {
        return  ln_3 + 2 * boost::math::digamma<RealType>(haplotype_pseudo_counts.at(genotype.at(0))) +
            boost::math::digamma<RealType>(haplotype_pseudo_counts.at((genotype.at(1) == genotype.at(0)) ?
                                                                      genotype.at(2) : genotype.at(1))) - k;
    } else if (genotype.num_occurences(genotype.at(1)) == 2) {
        return  ln_3 + 2 * boost::math::digamma<RealType>(haplotype_pseudo_counts.at(genotype.at(1))) +
            boost::math::digamma<RealType>(haplotype_pseudo_counts.at(genotype.at(0))) - k;
    } else {
        return ln_6 + boost::math::digamma<RealType>(haplotype_pseudo_counts.at(genotype.at(0))) +
            boost::math::digamma<RealType>(haplotype_pseudo_counts.at(genotype.at(1))) +
            boost::math::digamma<RealType>(haplotype_pseudo_counts.at(genotype.at(2))) + k;
    }
}

VariationalBayesGenotypeModel::RealType
VariationalBayesGenotypeModel::log_expected_genotype_probability_polyploid(const Genotype& genotype,
                                                                           const HaplotypePseudoCounts& haplotype_pseudo_counts) const
{
    //TODO
    return 0;
}

VariationalBayesGenotypeModel::HaplotypePseudoCounts
get_prior_pseudo_counts(const HaplotypePriors& the_haplotype_priors,
                        const Haplotype& the_reference_haplotype,
                        VariationalBayesGenotypeModel::RealType the_reference_haplotype_pseudo_count)
{
    VariationalBayesGenotypeModel::HaplotypePseudoCounts result {};
    result.reserve(the_haplotype_priors.size());
    
    static const VariationalBayesGenotypeModel::RealType ln_2 {std::log(2)};
    
    auto norm = sum(the_haplotype_priors);
    
    for (const auto& haplotype_prior : the_haplotype_priors) {
        result.emplace(haplotype_prior.first, digamma_inv(std::log(haplotype_prior.second / norm) -
                                                          std::log(the_reference_haplotype_pseudo_count - 0.5) -
                                                          ln_2));
    }
    
    result[the_reference_haplotype] = the_reference_haplotype_pseudo_count;
    
    return result;
}
