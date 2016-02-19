//
//  haplotype_phaser.cpp
//  Octopus
//
//  Created by Daniel Cooke on 30/04/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "haplotype_phaser.hpp"

#include <algorithm>
#include <iterator>
#include <cmath>
#include <functional>
#include <unordered_map>
#include <cassert>

#include "allele.hpp"
#include "mappable_algorithms.hpp"
#include "map_utils.hpp"

#include <iostream> // DEBUG

namespace Octopus
{

// public methods

auto calculate_max_indcluded(const unsigned max_haplotypes)
{
    return static_cast<unsigned>(std::max(1.0, std::log2(max_haplotypes)));
}

HaplotypePhaser::HaplotypePhaser(ContigNameType contig,
                                 const ReferenceGenome& reference,
                                 const std::vector<Variant>& candidates,
                                 const ReadMap& reads,
                                 const unsigned max_haplotypes,
                                 const unsigned max_indicators)
:
tree_ {contig, reference},
walker_ {max_indicators, calculate_max_indcluded(max_haplotypes)},
buffered_candidates_ {std::cbegin(candidates), std::cend(candidates)},
reads_ {&reads},
is_phasing_enabled_ {max_indicators != 0},
current_region_ {shift(head_region(buffered_candidates_.leftmost()), -1)},
next_region_ {walker_.walk(current_region_, *reads_, buffered_candidates_)},
max_haplotypes_ {max_haplotypes}
{}

bool HaplotypePhaser::done() const noexcept
{
    return buffered_candidates_.empty();
}

std::pair<std::vector<Haplotype>, GenomicRegion> HaplotypePhaser::fetch_next_haplotypes()
{
    std::swap(next_region_, current_region_);
    
    extend_tree(buffered_candidates_.overlap_range(current_region_), tree_);
    
    const auto haplotype_region = expand(current_region_, 50);
    
    return std::make_pair(tree_.extract_haplotypes(haplotype_region), current_region_);
}

void HaplotypePhaser::keep_haplotypes(const std::vector<Haplotype>& haplotypes)
{
    for (const auto& haplotype : haplotypes) {
        tree_.prune_unique(haplotype);
    }
}

void HaplotypePhaser::prepare_for_phasing(const HaplotypePosteiorMap& haplotype_posteriors)
{
    next_region_ = walker_.walk(current_region_, *reads_, buffered_candidates_);
    
    assert(current_region_ < next_region_);
    
    //std::cout << "next region is " << next_region_ << std::endl;
    
    if (!is_phasing_enabled() || is_empty_region(next_region_)) {
        tree_.clear();
        return;
    }
    
    const auto phased_region = left_overhang_region(current_region_, next_region_);
    
    tree_.remove(phased_region);
    
    if (tree_.empty()) return;
    
    using Pair = std::pair<HaplotypeReference, double>;
    
    std::vector<Pair> flat_haplotypes {std::cbegin(haplotype_posteriors), std::cend(haplotype_posteriors)};
    
    std::sort(std::begin(flat_haplotypes), std::end(flat_haplotypes),
              [] (const Pair& lhs, const Pair& rhs) {
                  return lhs.second < rhs.second;
              });
    
    const auto it = std::lower_bound(std::cbegin(flat_haplotypes), std::cend(flat_haplotypes),
                                     min_haplotype_posterior_,
                                     [] (const Pair& p, const double value) {
                                         return p.second < value;
                                     });
    
    //std::cout << "removing " << std::distance(std::cbegin(flat_haplotypes), it) << " haplotypes" << '\n';
    
    std::for_each(std::cbegin(flat_haplotypes), it,
                  [this] (const Pair& p) { tree_.prune_all(p.first); });
    
    const auto num_haplotypes_to_remove = calculate_num_haplotypes_to_remove();
    
    assert(num_haplotypes_to_remove < tree_.num_haplotypes());
    
    std::for_each(it, std::next(it, num_haplotypes_to_remove),
                  [this] (const Pair& p) { tree_.prune_all(p.first); });
}

boost::optional<HaplotypePhaser::PhaseSet>
HaplotypePhaser::phase(const std::vector<Haplotype>& haplotypes,
                       const GenotypePosteriors& genotype_posteriors)
{
    std::cout << "phasing" << std::endl;
    
    const auto phased_region = left_overhang_region(current_region_, next_region_);
    
    if (is_empty_region(phased_region)) return boost::none;
    
    auto variants = copy_contained(buffered_candidates_, phased_region);
    
    buffered_candidates_.erase_contained(phased_region);
    
    return find_optimal_phase_set(phased_region, std::move(variants), genotype_posteriors);
}

// private methods

unsigned HaplotypePhaser::calculate_num_haplotypes_to_remove() const
{
    const auto new_candidate_region = right_overhang_region(next_region_, current_region_);
    
    assert(!is_empty_region(new_candidate_region));
    
    const auto num_new_candidates = buffered_candidates_.count_overlapped(new_candidate_region);
    
    assert(num_new_candidates != 0);
    
    const auto cur_num_haplotypes = tree_.num_haplotypes();
    
    auto max_num_haplotypes_to_keep = max_haplotypes_ / std::exp2(num_new_candidates);
    
    if (cur_num_haplotypes <= max_num_haplotypes_to_keep) {
        return 0;
    } else {
        return cur_num_haplotypes - max_num_haplotypes_to_keep;
    }
}

bool HaplotypePhaser::is_phasing_enabled() const noexcept
{
    return is_phasing_enabled_;
}

// non-member methods

namespace
{
using GenotypeReference  = std::reference_wrapper<const Genotype<Haplotype>>;
using PhaseComplementSet = std::vector<GenotypeReference>;

struct PhaseComplementEqual
{
    PhaseComplementEqual(const std::vector<GenomicRegion>& variant_regions)
    : variant_regions_ {variant_regions} {}
    
    bool operator()(const GenotypeReference& lhs, const GenotypeReference& rhs) const
    {
        if (lhs.get().ploidy() == 1) return false;
        
        return std::all_of(std::cbegin(variant_regions_.get()), std::cend(variant_regions_.get()),
                           [&] (const auto& region) {
                               return are_equal_in_region<Haplotype>(lhs.get(), rhs.get(), region);
                           });
    }
private:
    std::reference_wrapper<const std::vector<GenomicRegion>> variant_regions_;
};

using PhaseComplementSetMap = std::unordered_map<GenotypeReference, PhaseComplementSet,
                                               std::hash<GenotypeReference>,
                                            PhaseComplementEqual>;

void add_to_phase_complement_set(const Genotype<Haplotype>& genotype,
                                 PhaseComplementSetMap& curr_result,
                                 const std::vector<GenomicRegion>& variant_regions)
{
    const auto it = curr_result.find(genotype);
    if (it == std::end(curr_result)) {
        curr_result[genotype].emplace_back(genotype);
    } else {
        it->second.emplace_back(genotype);
    }
}

using PhaseComplementSets = std::vector<PhaseComplementSet>;

PhaseComplementSets
partition_phase_complements(const std::vector<Genotype<Haplotype>>& genotypes,
                            const MappableSet<Variant>& variants)
{
    if (genotypes.empty()) return {};
    
    const auto variant_regions = extract_regions(variants);
    
    PhaseComplementSetMap result_map {genotypes.size(), std::hash<GenotypeReference> {},
                                        PhaseComplementEqual {variant_regions}};
    
    result_map[genotypes.front()].emplace_back(genotypes.front());
    
    std::for_each(std::next(std::cbegin(genotypes)), std::cend(genotypes),
                  [&] (const auto& genotype) {
                      add_to_phase_complement_set(genotype, result_map, variant_regions);
                  });
    
    PhaseComplementSets result {};
    result.reserve(result_map.size());
    
    for (auto&& p : result_map) {
        result.emplace_back(std::move(p.second));
    }
    
    return result;
}

double marginalise(const PhaseComplementSet& phase_set,
                   const HaplotypePhaser::SampleGenotypePosteriors& genotype_posteriors)
{
    return std::accumulate(std::cbegin(phase_set), std::cend(phase_set), 0.0,
                           [&] (const auto curr, const auto& genotype) {
                               return curr + genotype_posteriors.at(genotype);
                           });
}

double calculate_entropy(const PhaseComplementSet& phase_set,
                         const HaplotypePhaser::SampleGenotypePosteriors& genotype_posteriors)
{
    const auto norm = marginalise(phase_set, genotype_posteriors);
    return std::max(0.0, -std::accumulate(std::cbegin(phase_set), std::cend(phase_set), 0.0,
                                          [&genotype_posteriors, norm] (const auto curr, const auto& genotype) {
                                              const auto p = genotype_posteriors.at(genotype) / norm;
                                              return curr + p * std::log2(p);
                                          }));
}

double maximum_entropy(const size_t num_elements)
{
    return std::log2(num_elements);
}

double calculate_relative_entropy(const PhaseComplementSet& phase_set,
                                  const HaplotypePhaser::SampleGenotypePosteriors& genotype_posteriors)
{
    if (phase_set.size() < 2) return 1.0;
    return 1.0 - calculate_entropy(phase_set, genotype_posteriors) / maximum_entropy(phase_set.size());
}

double calculate_phase_score(const PhaseComplementSet& phase_set,
                             const HaplotypePhaser::SampleGenotypePosteriors& genotype_posteriors)
{
    return marginalise(phase_set, genotype_posteriors) * calculate_relative_entropy(phase_set, genotype_posteriors);
}

double calculate_phase_score(const PhaseComplementSets& phase_sets,
                             const HaplotypePhaser::SampleGenotypePosteriors& genotype_posteriors)
{
    return std::accumulate(std::cbegin(phase_sets), std::cend(phase_sets), 0.0,
                           [&] (auto curr, const auto& phase_set) {
                               return curr + calculate_phase_score(phase_set, genotype_posteriors);
                           });
}

HaplotypePhaser::SampleGenotypePosteriors
marginalise(const HaplotypePhaser::SampleGenotypePosteriors& genotype_posteriors, const GenomicRegion& region)
{
    HaplotypePhaser::SampleGenotypePosteriors result {};
    result.reserve(genotype_posteriors.size());
    
    for (const auto& p : genotype_posteriors) {
        result[splice<Haplotype>(p.first, region)] += p.second;
    }
    
    result.rehash(result.size());
    
    return result;
}

// TODO: this algorithm is very naive/slow right now.. needs much improvement
HaplotypePhaser::PhaseSet::SamplePhaseRegions
find_optimal_phase_regions(const GenomicRegion& region,
                           MappableSet<Variant> variants,
                           const HaplotypePhaser::SampleGenotypePosteriors& genotype_posteriors,
                           const double min_phase_score = 0.95)
{
    auto genotypes = extract_keys(genotype_posteriors);
    
    auto phase_set = partition_phase_complements(genotypes, variants);
    
    auto phase_score = calculate_phase_score(phase_set, genotype_posteriors);
    
    if (phase_score >= min_phase_score) {
        return {HaplotypePhaser::PhaseSet::PhaseRegion {region, phase_score}};
    }
    
    HaplotypePhaser::PhaseSet::SamplePhaseRegions result {};
    
    auto phase_begin_itr = std::cbegin(variants);
    auto phase_end_itr   = std::next(phase_begin_itr);
    
    double previous_phase_score {0};
    
    while (phase_begin_itr != std::cend(variants)) {
        MappableSet<Variant> curr_variants {phase_begin_itr, phase_end_itr};
        
        auto curr_region = encompassing_region(curr_variants);
        
        auto curr_genotype_posteriors = marginalise(genotype_posteriors, curr_region);
        
        genotypes = extract_keys(curr_genotype_posteriors);
        phase_set = partition_phase_complements(genotypes, curr_variants);
        
        phase_score = calculate_phase_score(phase_set, curr_genotype_posteriors);
        
        if (phase_score + 0.1 < previous_phase_score) {
            result.emplace_back(encompassing_region(phase_begin_itr, std::prev(phase_end_itr)),
                                previous_phase_score);
            phase_begin_itr = std::prev(phase_end_itr);
            previous_phase_score = 0.0;
        } else {
            if (phase_end_itr == std::cend(variants)) {
                result.emplace_back(curr_region, phase_score);
                phase_begin_itr = phase_end_itr;
            } else {
                previous_phase_score = phase_score;
                ++phase_end_itr;
            }
        }
    }
    
    return result;
}
} // namespace

HaplotypePhaser::PhaseSet
HaplotypePhaser::find_optimal_phase_set(const GenomicRegion& region,
                                        MappableSet<Variant> variants,
                                        const HaplotypePhaser::GenotypePosteriors& genotype_posteriors)
{
    HaplotypePhaser::PhaseSet result {region};
    
    if (variants.empty()) return result;
    
    result.phase_regions.reserve(genotype_posteriors.size());
    
    for (const auto& p : genotype_posteriors) {
        result.phase_regions.emplace(p.first, find_optimal_phase_regions(region, variants, p.second, min_phase_score_));
    }
    
    return result;
}

ProbabilityMatrix<GenomicRegion>
calculate_phase_scores(const ProbabilityMatrix<Genotype<Haplotype>>& genotype_posteriors)
{
    ProbabilityMatrix<GenomicRegion> result {};
    return result;
}
    
HaplotypePhaser::PhaseSet
phase(const GenomicRegion& region,
      const MappableSet<Variant>& variants,
      const ProbabilityMatrix<Genotype<Haplotype>>& genotype_posteriors)
{
    HaplotypePhaser::PhaseSet result {region};
    
    if (variants.empty()) return result;
    
    result.phase_regions.reserve(genotype_posteriors.size1());
    
    
    
    return result;
}
    
} // namespace Octopus
