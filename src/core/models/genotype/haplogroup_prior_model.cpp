// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "haplogroup_prior_model.hpp"

namespace octopus {

HaplogroupGenotypePriorModel::HaplogroupGenotypePriorModel(const GenotypePriorModel& germline_model,
                                                           SomaticMutationModel mutation_model)
: germline_model_ {germline_model}
, mutation_model_ {mutation_model}
{}

const GenotypePriorModel& HaplogroupGenotypePriorModel::germline_model() const noexcept
{
    return germline_model_.get();
}

SomaticMutationModel& HaplogroupGenotypePriorModel::mutation_model() noexcept
{
    return mutation_model_;
}

const SomaticMutationModel& HaplogroupGenotypePriorModel::mutation_model() const noexcept
{
    return mutation_model_;
}

HaplogroupGenotypePriorModel::HaplogroupPhylogeny
HaplogroupGenotypePriorModel::find_map_phylogeny(const Genotype<IndexedHaplotype<>>& genotype) const
{
    HaplogroupPhylogeny result {{genotype[0]}};
    for (unsigned k {1}; k < genotype.ploidy(); ++k) {
        result.add_descendant({genotype[k]}, genotype[0]);
    }
    return result;
}

} // namespace octopus
