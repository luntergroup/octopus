// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef haplogroup_prior_model_hpp
#define haplogroup_prior_model_hpp

#include <functional>
#include <cmath>
#include <iterator>
#include <algorithm>
#include <cassert>

#include "core/types/haplotype.hpp"
#include "core/types/indexed_haplotype.hpp"
#include "core/types/genotype.hpp"
#include "core/types/partitioned_genotype.hpp"
#include "core/types/phylogeny.hpp"
#include "genotype_prior_model.hpp"
#include "../mutation/somatic_mutation_model.hpp"
#include "utils/maths.hpp"

namespace octopus {

/*
The HaplogroupGenotypePriorModel puts a prior on a PartitionedGenotype with k partitions,
where each partition is a haplogroup. All haplotypes within a haplogroup have a common
acestor, that may or may not be observed in the haplogroup. Essentially, we're putting a prior
on a partially observed - incomplete - haploid phylogeny where a haplogroup is a clade on the tree. 
Since we won't be specifying the complete phylogeny within haplogroups, we (approximately) marginalise 
over all possible sub-phylogenies.
*/

class HaplogroupGenotypePriorModel
{
public:
    using LogProbability = double;
    
    HaplogroupGenotypePriorModel() = delete;
    
    HaplogroupGenotypePriorModel(const GenotypePriorModel& germline_model,
                                 SomaticMutationModel mutation_model);
    
    HaplogroupGenotypePriorModel(const HaplogroupGenotypePriorModel&)            = default;
    HaplogroupGenotypePriorModel& operator=(const HaplogroupGenotypePriorModel&) = default;
    HaplogroupGenotypePriorModel(HaplogroupGenotypePriorModel&&)                 = default;
    HaplogroupGenotypePriorModel& operator=(HaplogroupGenotypePriorModel&&)      = default;
    
    ~HaplogroupGenotypePriorModel() = default;
    
    const GenotypePriorModel& germline_model() const noexcept;
    SomaticMutationModel& mutation_model() noexcept;
    const SomaticMutationModel& mutation_model() const noexcept;
    
    template <std::size_t K>
    LogProbability evaluate(const PartitionedGenotype<IndexedHaplotype<>, K>& genotype) const;

private:
    std::reference_wrapper<const GenotypePriorModel> germline_model_;
    SomaticMutationModel mutation_model_;

    using HaplogroupPhylogeny = Phylogeny<IndexedHaplotype<>>;

    HaplogroupPhylogeny find_map_phylogeny(const Genotype<IndexedHaplotype<>>& genotype) const;
};

template <std::size_t K>
HaplogroupGenotypePriorModel::LogProbability 
HaplogroupGenotypePriorModel::evaluate(const PartitionedGenotype<IndexedHaplotype<>, K>& genotype) const
{
    LogProbability result {};
    Genotype<IndexedHaplotype<>> roots {K};
    for (std::size_t k {0}; k < K; ++k) {
        // Approximate full probability with most probable - a priori - phylogeny.
        // This works because we're assuming that the likelihood is independent of the
        // sub-phylogeny.
        const auto phylogeny = find_map_phylogeny(genotype[k]);
        for (const auto& haplotype : genotype.partition(k)) {
            if (haplotype == phylogeny.founder().id) {
                roots.emplace(haplotype);
            } else {
                result += mutation_model().evaluate(phylogeny.ancestor(haplotype).id, haplotype);
            }
        }
    }
    result += germline_model().evaluate(roots);
    return result;
}

} // namespace octopus

#endif
