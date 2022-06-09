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
    
    LogProbability evaluate(const PartitionedGenotype<IndexedHaplotype<>>& genotype) const;

private:
    std::reference_wrapper<const GenotypePriorModel> germline_model_;
    SomaticMutationModel mutation_model_;

    using HaplogroupPhylogeny = Phylogeny<IndexedHaplotype<>>;

    HaplogroupPhylogeny find_map_phylogeny(const Genotype<IndexedHaplotype<>>& genotype) const;
};

// non-member methods

template <typename Container1, typename Container2>
Container2&
evaluate(const Container1& genotypes, const HaplogroupGenotypePriorModel& model, Container2& result,
         const bool normalise = false, const bool add = false)
{
    if (add) {
        assert(result.size() == genotypes.size());
        std::transform(std::cbegin(genotypes), std::cend(genotypes), std::cbegin(result), std::begin(result),
                       [&] (const auto& genotype, auto curr) { return curr + model.evaluate(genotype); });
    } else {
        result.resize(genotypes.size());
        std::transform(std::cbegin(genotypes), std::cend(genotypes), std::begin(result),
                       [&] (const auto& genotype) { return model.evaluate(genotype); });
    }
    if (normalise) maths::normalise_logs(result);
    return result;
}

template <typename Container>
auto evaluate(const Container& genotypes, const HaplogroupGenotypePriorModel& model, const bool normalise = false)
{
    std::vector<HaplogroupGenotypePriorModel::LogProbability> result(genotypes.size());
    evaluate(genotypes, model, result, normalise, false);
    return result;
}

} // namespace octopus

#endif
