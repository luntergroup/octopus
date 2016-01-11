//
//  haplotype_phaser.cpp
//  Octopus
//
//  Created by Daniel Cooke on 30/04/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "haplotype_phaser.hpp"

#include <deque>
#include <algorithm>
#include <iterator>
#include <cmath>
#include <functional>

#include "allele.hpp"
#include "mappable_algorithms.hpp"
#include "map_utils.hpp"

#include <iostream> // DEBUG

namespace Octopus
{

// public methods

unsigned calculate_max_indcluded(unsigned max_haplotypes)
{
    return static_cast<unsigned>(std::max(1.0, std::log2(max_haplotypes)));
}

HaplotypePhaser::HaplotypePhaser(const ReferenceGenome& reference, const std::vector<Variant>& candidates,
                                 const ReadMap& reads, unsigned max_haplotypes, unsigned max_indicators)
:
tree_ {reference},
buffered_candidates_ {std::cbegin(candidates), std::cend(candidates)},
reads_ {&reads},
walker_ {max_indicators, calculate_max_indcluded(max_haplotypes),
    GenomeWalker::IndicatorLimit::SharedWithPreviousRegion, GenomeWalker::ExtensionLimit::SharedWithFrontier},
is_phasing_enabled_ {max_indicators != 0},
current_region_ {shift(get_head(buffered_candidates_.leftmost()), -1)},
next_region_ {walker_.walk(current_region_, *reads_, buffered_candidates_)}
{}

bool HaplotypePhaser::done() const noexcept
{
    return buffered_candidates_.empty();
}

std::vector<Haplotype> HaplotypePhaser::get_haplotypes()
{
    const auto next_candidates = buffered_candidates_.overlap_range(next_region_);
    
    for (const auto candidate : buffered_candidates_.overlap_range(next_region_)) {
        tree_.extend(candidate.get_ref_allele());
        tree_.extend(candidate.get_alt_allele());
    }
    
    current_region_ = next_region_;
    
    return tree_.get_haplotypes(current_region_);
}

std::vector<Haplotype> HaplotypePhaser::get_haplotypes(const GenotypePosteriors& genotype_posteriors)
{
    const auto next_candidates = buffered_candidates_.overlap_range(next_region_);
    
    for (const auto candidate : next_candidates) {
        tree_.extend(candidate.get_ref_allele());
        tree_.extend(candidate.get_alt_allele());
    }
    
    current_region_ = next_region_;
    
    return tree_.get_haplotypes(current_region_);
}

void HaplotypePhaser::set_haplotypes(const std::vector<Haplotype>& haplotypes)
{
    for (const auto& haplotype : haplotypes) {
        tree_.prune_unique(haplotype);
    }
}

boost::optional<HaplotypePhaser::PhaseSet>
HaplotypePhaser::phase(const std::vector<Haplotype>& haplotypes,
                       const GenotypePosteriors& genotype_posteriors)
{
    std::cout << "phasing" << std::endl;
    
    remove_low_posterior_haplotypes(haplotypes, genotype_posteriors);
    
    next_region_ = walker_.walk(current_region_, *reads_, buffered_candidates_);
    
    std::cout << "next region is " << next_region_ << std::endl;
    
    const auto phased_region = get_left_overhang(current_region_, next_region_);
    
    if (empty(phased_region)) {
        return boost::none;
    }
    
    std::cout << "phased region is " << phased_region << std::endl;
    
    auto variants = copy_overlapped(buffered_candidates_, phased_region);
    
    buffered_candidates_.erase_overlapped(phased_region);
    
    std::cout << "there are " << buffered_candidates_.size() << " candidates remaining" << std::endl;
    
    std::cout << "there are " << tree_.num_haplotypes() << " haplotypes in the tree" << std::endl;
    
    std::cout << "clearing " << phased_region << " from tree" << std::endl;
    
    tree_.remove(phased_region);
    
    if (tree_.empty()) {
        std::cout << "tree is empty" << std::endl;
    } else {
        std::cout << "there are " << tree_.num_haplotypes() << " haplotypes in the tree" << std::endl;
    }
    
    std::cout << "finding phase sets" << std::endl;
    
    auto result = find_optimal_phase_set(phased_region, std::move(variants), genotype_posteriors);
    
    std::cout << "done finding phase sets" << std::endl;
    
    return result;
}

// private methods

    std::unordered_map<std::reference_wrapper<const Haplotype>, double>
    compute_haplotype_posteriors(const std::vector<Haplotype>& haplotypes,
                                 const HaplotypePhaser::SampleGenotypePosteriors& genotype_posteriors)
    {
        std::unordered_map<std::reference_wrapper<const Haplotype>, double> result {};
        result.reserve(haplotypes.size());
        
        for (const auto& haplotype : haplotypes) {
            for (const auto& genotype_posterior: genotype_posteriors) {
                if (genotype_posterior.first.contains(haplotype)) {
                    result[haplotype] += genotype_posterior.second;
                }
            }
        }
        
        return result;
    }
    
    void HaplotypePhaser::remove_low_posterior_haplotypes(const std::vector<Haplotype>& haplotypes,
                                                          const GenotypePosteriors& genotype_posteriors)
    {
        if (!is_phasing_enabled_) {
            tree_.clear();
            return;
        }
        
        std::deque<std::reference_wrapper<const Haplotype>> low_posterior_haplotypes {};
        
        for (const auto& sample_genotype_posteriors : genotype_posteriors) {
            auto sample_haplotype_posteriors = compute_haplotype_posteriors(haplotypes, sample_genotype_posteriors.second);
            
            for (const auto& haplotype_posterior : sample_haplotype_posteriors) {
                if (haplotype_posterior.second < 1e-4) {
                    //print_variant_alleles(haplotype_posterior.first);
                    //std::cout << " " << haplotype_posterior.second << std::endl;
                    low_posterior_haplotypes.push_back(haplotype_posterior.first);
                }
            }
        }
        
        std::sort(std::begin(low_posterior_haplotypes), std::end(low_posterior_haplotypes));
        
        low_posterior_haplotypes.erase(std::unique(std::begin(low_posterior_haplotypes),
                                                   std::end(low_posterior_haplotypes)),
                                       std::end(low_posterior_haplotypes));
        
        std::cout << "removing " << low_posterior_haplotypes.size() << " haplotypes" << std::endl;
        
        for (const auto& haplotype : low_posterior_haplotypes) {
            tree_.prune_all(haplotype);
        }
    }
    
// non-member methods
    
    using PhaseComplementSet  = std::vector<Genotype<Haplotype>>;
    using PhaseComplementSets = std::vector<PhaseComplementSet>;
    
    PhaseComplementSets partition_phase_complements(const std::vector<Genotype<Haplotype>>& genotypes,
                                                    const std::vector<Variant>& variants);
    
    bool are_phase_complements(const Genotype<Haplotype>& lhs, const Genotype<Haplotype>& rhs,
                               const MappableSet<Variant>& variants)
    {
        if (lhs.ploidy() == 1) return false;
        
        const auto regions = get_regions(variants);
        
        return std::all_of(std::cbegin(regions), std::cend(regions),
                           [&lhs, &rhs] (const auto& region) {
                               return are_equal_in_region(lhs, rhs, region);
                           });
    }
    
    PhaseComplementSets
    partition_phase_complements(const std::vector<Genotype<Haplotype>>& genotypes,
                                const MappableSet<Variant>& variants)
    {
        using std::cbegin; using std::cend; using std::begin; using std::end; using std::for_each;
        
        PhaseComplementSets result {};
        
        if (genotypes.empty()) return result;
        
        result.reserve(genotypes.size());
        
        result.emplace_back(1, genotypes.front());
        
        for_each(std::next(cbegin(genotypes)), cend(genotypes),
                 [&result, &variants] (const auto& genotype) {
                     auto it = std::find_if(begin(result), end(result),
                                            [&genotype, &variants] (const auto& complements) {
                                                return are_phase_complements(genotype, complements.front(), variants);
                                            });
                     
                     if (it == end(result)) {
                         result.emplace_back(1, genotype);
                     } else {
                         it->emplace_back(genotype);
                     }
                 });
        
        return result;
    }
    
    double marginalise(const PhaseComplementSet& phase_set,
                       const HaplotypePhaser::SampleGenotypePosteriors& genotype_posteriors)
    {
        return std::accumulate(std::cbegin(phase_set), std::cend(phase_set), 0.0,
                               [&genotype_posteriors] (auto curr, const auto& genotype) {
                                   return curr + genotype_posteriors.at(genotype);
                               });
    }
    
    double calculate_entropy(const PhaseComplementSet& phase_set,
                             const HaplotypePhaser::SampleGenotypePosteriors& genotype_posteriors)
    {
        const auto norm = marginalise(phase_set, genotype_posteriors);
        return std::max(0.0, -std::accumulate(std::cbegin(phase_set), std::cend(phase_set), 0.0,
                                              [&genotype_posteriors, norm] (auto curr, const auto& genotype) {
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
                               [&genotype_posteriors] (auto curr, const auto& phase_set) {
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
        
        result.reserve(result.size());
        
        return result;
    }
    
    // TODO: this algorithm is very naive/slow right now.. needs much improvement
    HaplotypePhaser::PhaseSet::SamplePhaseRegions
    find_optimal_phase_regions(const GenomicRegion& region,
                               MappableSet<Variant> variants,
                               const HaplotypePhaser::SampleGenotypePosteriors& genotype_posteriors)
    {
        auto genotypes = get_keys(genotype_posteriors);
        
        auto phase_set = partition_phase_complements(genotypes, variants);
        
        auto phase_score = calculate_phase_score(phase_set, genotype_posteriors);
        
        std::cout << "phase score is " << phase_score << std::endl;
        
        if (phase_score >= 0.95) {
            return {HaplotypePhaser::PhaseSet::PhaseRegion {region, phase_score}};
        }
        
        HaplotypePhaser::PhaseSet::SamplePhaseRegions result {};
        
        auto phase_begin_itr = std::cbegin(variants);
        auto phase_end_itr   = std::next(phase_begin_itr);
        
        double previous_phase_score {};
        
        while (phase_begin_itr != std::cend(variants)) {
            MappableSet<Variant> curr_variants {phase_begin_itr, phase_end_itr};
            
            auto curr_region = get_encompassing_region(curr_variants);
            
            auto curr_genotype_posteriors = marginalise(genotype_posteriors, curr_region);
            
            genotypes = get_keys(curr_genotype_posteriors);
            phase_set = partition_phase_complements(genotypes, curr_variants);
            
            phase_score = calculate_phase_score(phase_set, curr_genotype_posteriors);
            
            std::cout << "current phase region is " << curr_region << std::endl;
            std::cout << "current phase score is " << phase_score << std::endl;
            
            if (phase_score + 0.1 < previous_phase_score) {
                auto phase_region = get_encompassing_region(phase_begin_itr, std::prev(phase_end_itr));
                std::cout << "adding phase region " << phase_region << std::endl;
                result.emplace_back(phase_region, previous_phase_score);
                phase_begin_itr = std::prev(phase_end_itr);
                previous_phase_score = 0.0;
            } else {
                if (phase_end_itr == std::cend(variants)) {
                    std::cout << "adding last phase region" << std::endl;
                    result.emplace_back(curr_region, phase_score);
                    phase_begin_itr = phase_end_itr;
                } else {
                    std::cout << "moving on" << std::endl;
                    previous_phase_score = phase_score;
                    ++phase_end_itr;
                }
            }
        }
        
        return result;
    }
    
    HaplotypePhaser::PhaseSet
    HaplotypePhaser::find_optimal_phase_set(const GenomicRegion& region,
                                            MappableSet<Variant> variants,
                                            const HaplotypePhaser::GenotypePosteriors& genotype_posteriors)
    {
        HaplotypePhaser::PhaseSet result {region};
        
        if (variants.empty()) return result;
        
        result.phase_regions.reserve(genotype_posteriors.size());
        
        for (const auto& p : genotype_posteriors) {
            result.phase_regions.emplace(p.first, find_optimal_phase_regions(region, variants, p.second));
        }
        
        return result;
    }

} // namespace Octopus
