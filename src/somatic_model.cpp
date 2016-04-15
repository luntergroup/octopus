//
//  somatic_model.cpp
//  Octopus
//
//  Created by Daniel Cooke on 12/04/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#include "somatic_model.hpp"

#include <utility>
#include <cmath>
#include <numeric>

namespace Octopus
{
    SomaticModel::SomaticModel(const CoalescentModel& germline_model, double somatic_mutation_rate)
    :
    germline_model_ {germline_model},
    somatic_mutation_rate_ {somatic_mutation_rate}
    {}
    
    const CoalescentModel& SomaticModel::get_germline_model() const noexcept
    {
        return germline_model_;
    }
    
    double probability_of_somatic(const Haplotype& somatic, const Haplotype& germline,
                                  double somatic_mutation_probability = 0.00001)
    {
        const auto variants = difference(somatic, germline);
        return std::pow(somatic_mutation_probability, variants.size());
    }
    
    // p(somatic | germline) = 1 / M sum k = 1 -> M p(somatic | germline_k) (M = germline ploidy)
    double probability_of_somatic(const Haplotype& somatic,
                                  const Genotype<Haplotype>& germline_genotype,
                                  const CoalescentModel& germline_model)
    {
        return std::accumulate(std::cbegin(germline_genotype), std::cend(germline_genotype), 0.0,
                               [&somatic] (const auto curr, const Haplotype& germline) {
                                   return curr + probability_of_somatic(somatic, germline);
                               }) / germline_genotype.ploidy();
    }
    
    double SomaticModel::evaluate(const CancerGenotype<Haplotype>& genotype) const
    {
        const auto& germline = genotype.get_germline_genotype();
        const auto& somatic  = genotype.get_cancer_element();
        
        const auto germline_log_prior = germline_model_.get().evaluate(germline);
        
        const auto somatic_probability_given_germline = probability_of_somatic(somatic, germline,
                                                                               germline_model_);
        
        return germline_log_prior + std::log(somatic_probability_given_germline);
    }
} // namespace Octopus
