// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "cancer_genotype_prior_model.hpp"

#include <utility>
#include <numeric>
#include <stdexcept>

#include "utils/mappable_algorithms.hpp"

namespace octopus {

CancerGenotypePriorModel::CancerGenotypePriorModel(const GenotypePriorModel& germline_model,
                                                   SomaticMutationModel mutation_model)
: germline_model_ {germline_model}
, mutation_model_ {mutation_model}
{}

const GenotypePriorModel& CancerGenotypePriorModel::germline_model() const noexcept
{
    return germline_model_.get();
}

SomaticMutationModel& CancerGenotypePriorModel::mutation_model() noexcept
{
    return mutation_model_;
}

const SomaticMutationModel& CancerGenotypePriorModel::mutation_model() const noexcept
{
    return mutation_model_;
}

double CancerGenotypePriorModel::evaluate(const CancerGenotype<Haplotype>& genotype) const
{
    const auto& germline = genotype.germline_genotype();
    const auto& somatic = genotype.somatic_element();
    const auto germline_log_prior = germline_model_.get().evaluate(germline);
    return germline_log_prior + ln_probability_of_somatic_given_genotype(somatic, germline);
}

double CancerGenotypePriorModel::evaluate(const std::vector<unsigned>& germline_indices, const unsigned somatic_index) const
{
    return germline_model_.get().evaluate(germline_indices)
           + ln_probability_of_somatic_given_genotype(somatic_index, germline_indices);
}

double CancerGenotypePriorModel::ln_probability_of_somatic_given_haplotype(const Haplotype& somatic, const Haplotype& germline) const
{
    return mutation_model_.evaluate(somatic, germline);
}

double CancerGenotypePriorModel::ln_probability_of_somatic_given_haplotype(unsigned somatic_index, unsigned germline_index) const
{
    return mutation_model_.evaluate(somatic_index, germline_index);
}

} // namespace octopus
