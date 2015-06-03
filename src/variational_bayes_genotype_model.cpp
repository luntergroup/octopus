//
//  variational_bayes_genotype_model.cpp
//  Octopus
//
//  Created by Daniel Cooke on 01/04/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "variational_bayes_genotype_model.h"

#include <boost/math/special_functions/digamma.hpp>

namespace Octopus
{

namespace BayesianGenotypeModel
{
    VariationalBayesGenotypeModel::VariationalBayesGenotypeModel(ReadModel& read_model, unsigned ploidy)
    :
    ploidy_ {ploidy},
    read_model_ {read_model}
    {}
    
    VariationalBayesGenotypeModel::RealType
    VariationalBayesGenotypeModel::log_expected_genotype_probability(const Genotype<Haplotype>& genotype,
                                                                     const HaplotypePseudoCounts<RealType>& haplotype_pseudo_counts)
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
                                                            const SampleGenotypeProbabilities<RealType>& sample_genotype_responsabilities) const
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
                                                                    const GenotypeProbabilities<SampleIdType, RealType>& genotype_responsabilities) const
    {
        RealType result {prior_pseudo_count};
        
        for (const auto& sample_genotype_responsabilities : genotype_responsabilities) {
            result += expected_haplotype_count(haplotype, sample_genotype_responsabilities.second);
        }
        
        return result;
    }
    
    void VariationalBayesGenotypeModel::clear_cache()
    {
        read_model_.clear_cache();
    }
    
    // Private methods
    
    VariationalBayesGenotypeModel::RealType
    VariationalBayesGenotypeModel::log_expected_genotype_probability_haploid(const Genotype<Haplotype>& genotype,
                                                                             const HaplotypePseudoCounts<RealType>& haplotype_pseudo_counts) const
    {
        return boost::math::digamma<RealType>(haplotype_pseudo_counts.at(genotype.at(0)))
                - boost::math::digamma<RealType>(sum_values(haplotype_pseudo_counts));
    }
    
    VariationalBayesGenotypeModel::RealType
    VariationalBayesGenotypeModel::log_expected_genotype_probability_diploid(const Genotype<Haplotype>& genotype,
                                                                             const HaplotypePseudoCounts<RealType>& haplotype_pseudo_counts) const
    {
        const static RealType ln_2 {std::log(2)};
        
        return ((genotype.is_homozygous()) ?
                2 * boost::math::digamma<RealType>(haplotype_pseudo_counts.at(genotype.at(0)))
                :
                ln_2 + boost::math::digamma<RealType>(haplotype_pseudo_counts.at(genotype.at(0))) +
                    boost::math::digamma<RealType>(haplotype_pseudo_counts.at(genotype.at(1))))
                    - 2 * boost::math::digamma<RealType>(sum_values(haplotype_pseudo_counts));
    }
    
    VariationalBayesGenotypeModel::RealType
    VariationalBayesGenotypeModel::log_expected_genotype_probability_triploid(const Genotype<Haplotype>& genotype,
                                                                              const HaplotypePseudoCounts<RealType>& haplotype_pseudo_counts) const
    {
        const static RealType ln_3 {std::log(3)};
        const static RealType ln_6 {std::log(6)};
        
        const auto k = 3 * boost::math::digamma<RealType>(sum_values(haplotype_pseudo_counts));
        
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
                        boost::math::digamma<RealType>(haplotype_pseudo_counts.at(genotype.at(2))) - k;
        }
    }
    
    VariationalBayesGenotypeModel::RealType
    VariationalBayesGenotypeModel::log_expected_genotype_probability_polyploid(const Genotype<Haplotype>& genotype,
                                                                               const HaplotypePseudoCounts<RealType>& haplotype_pseudo_counts) const
    {
        const auto k = 3 * boost::math::digamma<RealType>(sum_values(haplotype_pseudo_counts));
        
        auto unique_haplotypes = genotype.get_unique();
        
        RealType r {0};
        
        std::vector<unsigned> occurences {};
        occurences.reserve(unique_haplotypes.size());
        
        for (const auto& haplotype : unique_haplotypes) {
            occurences.push_back(genotype.num_occurences(haplotype));
            r += genotype.num_occurences(haplotype) * haplotype_pseudo_counts.at(haplotype);
        }
        
        return log_multinomial_coefficient<RealType>(occurences.cbegin(), occurences.cend()) + r - k;
    }
    
} // end namespace BayesianGenotypeModel

} // end namespace Octopus
