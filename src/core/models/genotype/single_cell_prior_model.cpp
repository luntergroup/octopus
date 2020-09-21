// Copyright (c) 2015-2020 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "single_cell_prior_model.hpp"

#include <utility>
#include <algorithm>
#include <numeric>
#include <iterator>
#include <cmath>

#include "utils/maths.hpp"

namespace octopus { namespace model {

SingleCellPriorModel::SingleCellPriorModel(CellPhylogeny phylogeny,
                                           const GenotypePriorModel& germline_prior_model,
                                           const DeNovoModel& denovo_model,
                                           Parameters parameters)
: phylogeny_ {std::move(phylogeny)}
, germline_prior_model_ {germline_prior_model}
, denovo_model_ {denovo_model}
, parameters_ {std::move(parameters)}
{}

const SingleCellPriorModel::CellPhylogeny& SingleCellPriorModel::phylogeny() const noexcept
{
    return phylogeny_;
}

const GenotypePriorModel& SingleCellPriorModel::germline_prior_model() const noexcept
{
    return germline_prior_model_.get();
}

const DeNovoModel& SingleCellPriorModel::denovo_model() const noexcept
{
    return denovo_model_.get();
}

SingleCellPriorModel::LogProbability
SingleCellPriorModel::evaluate(const std::vector<Genotype<IndexedHaplotype<>>>& genotypes) const
{
    LogProbability result {0};
    for (std::size_t id {0}; id < genotypes.size(); ++id) {
        if (id == phylogeny_.founder().id) {
            result += germline_prior_model_.get().evaluate(genotypes[id]);
        } else {
            result += log_probability(genotypes[phylogeny_.ancestor(id).id], genotypes[id]);
        }
    }
    return result;
}

SingleCellPriorModel::LogProbability
SingleCellPriorModel::evaluate(const std::vector<GenotypeReference>& genotypes) const
{
    LogProbability result {0};
    for (std::size_t id {0}; id < genotypes.size(); ++id) {
        if (id == phylogeny_.founder().id) {
            result += germline_prior_model_.get().evaluate(genotypes[id]);
        } else {
            result += log_probability(genotypes[phylogeny_.ancestor(id).id], genotypes[id]);
        }
    }
    return result;
}

SingleCellPriorModel::LogProbability
SingleCellPriorModel::compute_cnv_log_prior(const Genotype<IndexedHaplotype<>>& ancestor) const
{
    // This is the prior probability that an entire haplotype is gained or lost from the ancestor genotype
    // We assume that the length distribution is such that any copy number event will affect the whole haplotype.
    // Maybe for very large haplotypes we'll need to look at another model
    return std::log(parameters_.copy_number_prior);
}

SingleCellPriorModel::LogProbability
SingleCellPriorModel::log_probability(const Genotype<IndexedHaplotype<>>& ancestor, const Genotype<IndexedHaplotype<>>& descendant) const
{
    const auto cnv_log_prior = compute_cnv_log_prior(ancestor);
    LogProbability result {0};
    std::vector<unsigned> ancestor_haplotype_indices(ancestor.ploidy());
    std::iota(std::begin(ancestor_haplotype_indices), std::end(ancestor_haplotype_indices), 0u);
    if (ancestor.ploidy() != descendant.ploidy()) {
        const auto p = std::minmax({ancestor.ploidy(), descendant.ploidy()});
        const auto copy_change = p.second - p.first;
        result += cnv_log_prior * copy_change;
        if (ancestor.ploidy() < descendant.ploidy()) {
            // Copy gain
            for (unsigned i {0}; i < ancestor.ploidy(); ++i) {
                ancestor_haplotype_indices.insert(std::cend(ancestor_haplotype_indices), copy_change, i);
            }
        }
    }
    std::vector<LogProbability> conditional_log_probabilities {};
    conditional_log_probabilities.reserve(maths::factorial<std::size_t>(ancestor_haplotype_indices.size()));
    do {
        conditional_log_probabilities.push_back(
            std::inner_product(std::cbegin(descendant), std::cend(descendant), std::cbegin(ancestor_haplotype_indices),
                           0.0, std::plus<> {}, [&] (const auto& descendant_haplotype, const auto ancestor_haplotype_idx) {
            return denovo_model_.get().evaluate(descendant_haplotype, ancestor[ancestor_haplotype_idx]); }));
    } while (std::next_permutation(std::begin(ancestor_haplotype_indices), std::end(ancestor_haplotype_indices)));
    result += maths::log_sum_exp(conditional_log_probabilities) - std::log(conditional_log_probabilities.size());
    return result;
}

} // namespace model
} // namespace octopus
