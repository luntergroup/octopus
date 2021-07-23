// Copyright (c) 2015-2021 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "pedigree_prior_model.hpp"

#include <algorithm>
#include <iterator>
#include <utility>

namespace octopus {

namespace {

template <typename T>
std::size_t index_of(const T& value, const std::vector<T>& values)
{
    auto itr = std::find(std::cbegin(values), std::cend(values), value);
    return std::distance(std::cbegin(values), itr);
}

} // namespace

PedigreePriorModel::PedigreePriorModel(const SampleVector& samples,
                                       const Pedigree& pedigree,
                                       std::unique_ptr<PopulationPriorModel> population_prior_model,
                                       InheritanceModel inheritance_model)
: population_prior_model_ {std::move(population_prior_model)}
, inheritance_model_ {std::move(inheritance_model)}
{
    // TODO: add twins
    founders_.reserve(samples.size());
    offspring_with_two_parents_.reserve(samples.size());
    offspring_with_one_parent_.reserve(samples.size());
    for (std::size_t s {0}; s < samples.size(); ++s) {
        if (pedigree.is_founder(samples[s])) {
            founders_.push_back(s);
        } else if (pedigree.num_parents(samples[s]) == 2) {
            auto mother_idx = index_of(*pedigree.mother_of(samples[s]), samples);
            auto father_idx = index_of(*pedigree.father_of(samples[s]), samples);
            offspring_with_two_parents_.emplace_back(s, std::make_pair(mother_idx, father_idx));
        } else {
            auto parent = pedigree.mother_of(samples[s]);
            if (!parent) parent = pedigree.father_of(samples[s]);
            assert(parent);
            offspring_with_one_parent_.emplace_back(s, index_of(*parent, samples));
        }
    }
}

const PopulationPriorModel& PedigreePriorModel::population_model() const noexcept
{
    return *population_prior_model_;
}

const InheritanceModel& PedigreePriorModel::inheritance_model() const noexcept
{
    return inheritance_model_;
}

void PedigreePriorModel::prime(const MappableBlock<Haplotype>& haplotypes)
{
    population_prior_model_->prime(haplotypes);
    inheritance_model_.prime(haplotypes);
}

void PedigreePriorModel::unprime() noexcept
{
    population_prior_model_->unprime();
    inheritance_model_.unprime();
}

bool PedigreePriorModel::is_primed() const noexcept
{
    return population_prior_model_->is_primed() && inheritance_model_.is_primed();
}

PedigreePriorModel::LogProbability 
PedigreePriorModel::evaluate(const std::vector<Genotype<IndexedHaplotype<>>>& genotypes) const
{
    std::vector<Genotype<IndexedHaplotype<>>> population_genotypes {};
    population_genotypes.reserve(genotypes.size());
    for (auto s : founders_) {
        population_genotypes.push_back(genotypes[s]);
    }
    // For offspring without both parents, there are (for autosomes) haplotypes
    // not explained by any ancestors. These must be accounted for, and can be 
    // considered part of the background population haplotypes. The problem is that
    // we don't know which haplotypes are derived and which are background a-priori,
    // and evaluating all possabilities is not practical. We approximate the full
    // distribution by adding all single-parent offspring genotypes to the population
    // set. This mostly has the affect of penalising de novo mutations twice (once in
    // the population model and once in the de novo model).
    for (const auto& p : offspring_with_one_parent_) {
        population_genotypes.push_back(genotypes[p.first]);
    }
    for (const auto& p : twins_with_one_parent_) {
        population_genotypes.push_back(genotypes[p.first.first]);
        population_genotypes.push_back(genotypes[p.first.second]);
    }
    auto result = population_model().evaluate(population_genotypes);
    for (const auto& p : offspring_with_two_parents_) {
        result += inheritance_model().evaluate(genotypes[p.first], genotypes[p.second.first], genotypes[p.second.second]);
    }
    for (const auto& p : offspring_with_one_parent_) {
        result += inheritance_model().evaluate(genotypes[p.first], genotypes[p.second]);
    }
    for (const auto& p : twins_with_two_parents_) {
        result += inheritance_model().evaluate_twins(genotypes[p.first.first], genotypes[p.first.second], genotypes[p.second.first], genotypes[p.second.second]);
    }
    for (const auto& p : twins_with_one_parent_) {
        result += inheritance_model().evaluate_twins(genotypes[p.first.first], genotypes[p.first.second], genotypes[p.second]);
    }
    return result;
}

} // namespace octopus
