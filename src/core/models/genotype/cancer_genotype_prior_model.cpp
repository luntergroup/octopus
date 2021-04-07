// Copyright (c) 2015-2021 Daniel Cooke
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

CancerGenotypePriorModel::LogProbability CancerGenotypePriorModel::evaluate(const CancerGenotype<IndexedHaplotype<>>& genotype) const
{
    const auto& germline_genotype = genotype.germline();
    auto result = germline_model_.get().evaluate(germline_genotype);
    // Model assumes independence between somatic haplotypes given germline genotype
    for (const auto& somatic_haplotype : genotype.somatic()) {
        result += ln_probability_of_somatic_given_genotype(somatic_haplotype, germline_genotype);
    }
    return result;
}

CancerGenotypePriorModel::LogProbability
CancerGenotypePriorModel::ln_probability_of_somatic_given_haplotype(const IndexedHaplotype<>& somatic, const IndexedHaplotype<>& germline) const
{
    return mutation_model_.evaluate(somatic, germline);
}

} // namespace octopus
