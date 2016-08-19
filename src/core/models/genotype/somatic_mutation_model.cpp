// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "somatic_mutation_model.hpp"

#include <utility>
#include <cmath>
#include <numeric>
#include <stdexcept>

namespace octopus {

SomaticMutationModel::SomaticMutationModel(const CoalescentModel& germline_model,
                                           double somatic_mutation_rate)
: germline_model_ {germline_model}
, somatic_mutation_rate_ {somatic_mutation_rate}
{
    if (somatic_mutation_rate <= 0) {
        throw std::domain_error {"SomaticMutationModel: somatic mutation rate must be > 0"};
    }
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

double SomaticMutationModel::evaluate(const CancerGenotype<Haplotype>& genotype) const
{
    const auto& germline = genotype.germline_genotype();
    const auto& somatic  = genotype.somatic_element();
    
    const auto germline_log_prior = germline_model_.get().evaluate(germline);
    
    const auto somatic_probability_given_germline = probability_of_somatic(somatic, germline,
                                                                           germline_model_);
    
    return germline_log_prior + std::log(somatic_probability_given_germline);
}

} // namespace octopus
