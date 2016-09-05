// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "somatic_mutation_model.hpp"

#include <utility>
#include <cmath>
#include <numeric>
#include <stdexcept>

namespace octopus {

SomaticMutationModel::SomaticMutationModel(const CoalescentModel& germline_model,
                                           Parameters params)
: germline_model_ {germline_model}
, params_ {params}
{
    if (params_.somatic_mutation_rate <= 0) {
        throw std::domain_error {"SomaticMutationModel: somatic mutation rate must be > 0"};
    }
}

double SomaticMutationModel::evaluate(const CancerGenotype<Haplotype>& genotype) const
{
    const auto& germline = genotype.germline_genotype();
    const auto& somatic = genotype.somatic_element();
    const auto germline_log_prior = germline_model_.get().evaluate(germline);
    const auto somatic_probability_given_germline = probability_of_somatic(somatic, germline);
    return germline_log_prior + std::log(somatic_probability_given_germline);
}

// p(somatic | germline) = 1 / M sum k = 1 -> M p(somatic | germline_k) (M = germline ploidy)
double SomaticMutationModel::probability_of_somatic(const Haplotype& somatic, const Genotype<Haplotype>& germline) const
{
    return std::accumulate(std::cbegin(germline), std::cend(germline), 0.0,
                           [this, &somatic](const auto curr, const Haplotype& germline) {
                               return curr + probability_of_somatic(somatic, germline);
                           }) / germline.ploidy();
}

double SomaticMutationModel::probability_of_somatic(const Haplotype& somatic, const Haplotype& germline) const
{
    const auto variants = difference(somatic, germline);
    return std::pow(params_.somatic_mutation_rate, variants.size());
}

} // namespace octopus
