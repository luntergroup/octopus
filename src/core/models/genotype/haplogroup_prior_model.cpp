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

HaplogroupGenotypePriorModel::LogProbability 
HaplogroupGenotypePriorModel::evaluate(const PartitionedGenotype<IndexedHaplotype<>>& genotype) const
{
    LogProbability result {};
    const auto num_partions = static_cast<unsigned>(genotype.num_partitions());
    Genotype<IndexedHaplotype<>> roots {num_partions};
    for (unsigned k {0}; k < num_partions; ++k) {
        // Approximate full probability with most probable - a priori - phylogeny.
        // This works because we're assuming that the likelihood is independent of the
        // sub-phylogeny.
        const auto& partition = genotype.partition(k);
        if (partition.ploidy() > 0) {
            const auto phylogeny = find_map_phylogeny(partition);
            for (const auto& haplotype : partition) {
                if (haplotype == phylogeny.founder().id) {
                    roots.emplace(haplotype);
                } else {
                    result += mutation_model().evaluate(phylogeny.ancestor(haplotype).id, haplotype);
                }
            }
        }
    }
    result += germline_model().evaluate(roots);
    return result;
}

HaplogroupGenotypePriorModel::HaplogroupPhylogeny
HaplogroupGenotypePriorModel::find_map_phylogeny(const Genotype<IndexedHaplotype<>>& genotype) const
{
    HaplogroupPhylogeny result {{genotype[0]}};
    for (unsigned k {1}; k < genotype.ploidy(); ++k) {
        result.add_descendant({genotype[k]}, genotype[k - 1]);
    }
    return result;
}

} // namespace octopus
