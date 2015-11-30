//
//  haplotype_phaser.cpp
//  Octopus
//
//  Created by Daniel Cooke on 30/04/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "haplotype_phaser.hpp"

#include <deque>
#include <algorithm> // std::max, std::sort, std::unique, std::for_each
#include <iterator>  // std::begin, std::end, std::cbegin, std::cend
#include <cmath>     // std::log2

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

HaplotypePhaser::HaplotypePhaser(ReferenceGenome& reference, const std::vector<Variant>& candidates,
                                 const ReadMap& reads, unsigned max_haplotypes, unsigned max_indicators)
:
tree_ {reference},
buffered_candidates_ {std::cbegin(candidates), std::cend(candidates)},
reads_ {&reads},
walker_ {max_indicators, calculate_max_indcluded(max_haplotypes),
    GenomeWalker::IndicatorLimit::SharedWithPreviousRegion, GenomeWalker::ExtensionLimit::NoLimit},
is_phasing_enabled_ {max_indicators != 0},
tree_region_ {shift(get_head(buffered_candidates_.leftmost()), -1)},
next_region_ {walker_.walk(tree_region_, *reads_, buffered_candidates_)}
{}

bool HaplotypePhaser::done() const noexcept
{
    return buffered_candidates_.empty();
}

std::vector<Haplotype> HaplotypePhaser::get_haplotypes()
{
    auto contained_candidates = buffered_candidates_.contained_range(next_region_);
    
    std::for_each(std::cbegin(contained_candidates), std::cend(contained_candidates),
                  [this] (const auto& candidate) {
                      tree_.extend(candidate.get_reference_allele());
                      tree_.extend(candidate.get_alternative_allele());
                  });
    
    tree_region_ = next_region_;
    
    return tree_.get_haplotypes(tree_region_);
}

void HaplotypePhaser::unique(const std::vector<Haplotype>& haplotypes)
{
    for (const auto& haplotype : haplotypes) {
        tree_.prune_unique(haplotype);
    }
}

std::unordered_map<Haplotype, double>
compute_haplotype_posteriors(const std::vector<Haplotype>& haplotypes,
                             const HaplotypePhaser::SampleGenotypePosteriors& genotype_posteriors)
{
    std::unordered_map<Haplotype, double> result {};
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

HaplotypePhaser::PhaseSet
HaplotypePhaser::phase(const std::vector<Haplotype>& haplotypes,
                       const GenotypePosteriors& genotype_posteriors)
{
    next_region_ = walker_.walk(tree_region_, *reads_, buffered_candidates_);
    
    auto phased_region = get_left_overhang(tree_region_, next_region_);
    
    auto contained = buffered_candidates_.contained_range(phased_region);
    
    MappableSet<Variant> variants {contained.begin(), contained.end()};
    
    buffered_candidates_.erase_contained(phased_region);
    
    tree_.clear(phased_region);
    
    if (is_phasing_enabled_) {
        std::deque<Haplotype> low_posterior_haplotypes {};
        
        for (const auto& sample_genotype_posteriors : genotype_posteriors) {
            auto sample_haplotype_posteriors = compute_haplotype_posteriors(haplotypes, sample_genotype_posteriors.second);
            
            for (const auto& haplotype_posterior : sample_haplotype_posteriors) {
                if (haplotype_posterior.second < 0.0001) {
                    low_posterior_haplotypes.push_back(haplotype_posterior.first);
                }
            }
        }
        
        std::sort(std::begin(low_posterior_haplotypes), std::end(low_posterior_haplotypes));
        
        //std::cout << "removing " << low_posterior_haplotypes.size() << " haplotypes" << std::endl;
        
        low_posterior_haplotypes.erase(std::unique(std::begin(low_posterior_haplotypes),
                                                   std::end(low_posterior_haplotypes)),
                                       std::end(low_posterior_haplotypes));
        
        for (const auto& haplotype : low_posterior_haplotypes) {
            tree_.prune_all(haplotype);
        }
    } else {
        tree_.clear();
    }
    
    return find_optimal_phase_set(phased_region, std::move(variants), genotype_posteriors);
}

// private methods

// non-member methods
    
    using PhaseComplementSet  = std::vector<Genotype<Haplotype>>;
    using PhaseComplementSets = std::vector<PhaseComplementSet>;
    
    PhaseComplementSets partition_phase_complements(const std::vector<Genotype<Haplotype>>& genotypes,
                                                    const std::vector<Variant>& variants);
    
    bool are_phase_complements(const Genotype<Haplotype>& lhs, const Genotype<Haplotype>& rhs,
                               const MappableSet<Variant>& variants)
    {
        if (lhs.ploidy() == 1) return false;
        
        auto regions = get_regions(variants);
        
        return std::all_of(std::cbegin(regions), std::cend(regions),
                           [&lhs, &rhs] (const auto& region) {
                               return equal_in_region(lhs, rhs, region);
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
                                                return std::any_of(cbegin(complements), cend(complements),
                                                                   [&genotype, &variants] (const auto& complement) {
                                                                       return are_phase_complements(genotype, complement, variants);
                                                                   });
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
        return std::max(0.0, -std::accumulate(std::cbegin(phase_set), std::cend(phase_set), 0.0,
                                              [&genotype_posteriors] (auto curr, const auto& genotype) {
                                                  const auto p = genotype_posteriors.at(genotype);
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
    
    HaplotypePhaser::PhaseSet::SamplePhaseRegions
    find_optimal_phase_regions(const GenomicRegion& region,
                               MappableSet<Variant> variants,
                               const HaplotypePhaser::SampleGenotypePosteriors& genotype_posteriors)
    {
        auto genotypes = get_keys(genotype_posteriors);
        
        auto phase_set = partition_phase_complements(genotypes, variants);
        
        auto phase_score = calculate_phase_score(phase_set, genotype_posteriors);
        
        std::cout << "phase score is " << phase_score << std::endl;
        
        HaplotypePhaser::PhaseSet::PhaseRegion result {region, phase_score};
        
        return {result};
    }
    
    HaplotypePhaser::PhaseSet
    HaplotypePhaser::find_optimal_phase_set(const GenomicRegion& region,
                                            MappableSet<Variant> variants,
                                            const HaplotypePhaser::GenotypePosteriors& genotype_posteriors)
    {
        HaplotypePhaser::PhaseSet result {region};
        result.phase_regions.reserve(genotype_posteriors.size());
        
        for (const auto& p : genotype_posteriors) {
            result.phase_regions.emplace(p.first, find_optimal_phase_regions(region, variants, p.second));
        }
        
        return result;
    }

} // namespace Octopus
