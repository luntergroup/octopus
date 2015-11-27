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

#include <iostream> // DEBUG

namespace Octopus
{

// non-member methods

std::vector<std::vector<Genotype<Haplotype>>>
partition_phase_complements(const std::vector<Genotype<Haplotype>>& genotypes, const std::vector<Variant>& variants);

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
    
    std::vector<Variant> variants {contained.begin(), contained.end()};
    
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
    
    PhaseSet result {phased_region};
    
    for (const auto& p : genotype_posteriors) {
        const auto& sample_genotype_posteriors = p.second;
        
        std::vector<Genotype<Haplotype>> genotypes {};
        genotypes.reserve(genotype_posteriors.size());
        
        std::transform(std::cbegin(sample_genotype_posteriors), std::cend(sample_genotype_posteriors),
                       std::back_inserter(genotypes), [] (const auto& p) { return p.first; });
        
        auto phase_complements = partition_phase_complements(genotypes, variants);
        
        double phase_quality {1.0};
        
        for (const auto& set : phase_complements) {
            if (set.size() > 1) {
                auto c = std::count_if(std::cbegin(set), std::cend(set),
                                       [&sample_genotype_posteriors] (const auto& genotype) {
                                           return sample_genotype_posteriors.at(genotype) > 0.000001;
                                       });
                
                if (c > 1) {
                    std::cout << "phase complements set: " << std::endl;
                    for (const auto& g : set) {
                        print_variant_alleles(g);
                        std::cout << std::endl;
                    }
                    
                    auto entropy = -std::accumulate(std::cbegin(set), std::cend(set), 0.0,
                                                   [sample_genotype_posteriors] (auto curr, const auto& genotype) {
                                                       auto p = sample_genotype_posteriors.at(genotype);
                                                       return curr + (p * std::log2(p));
                                                   });
                    
                    phase_quality -= entropy;
                }
            }
        }
        
        if (phase_quality < 0.0) phase_quality = 0.0;
        
        std::cout << "phase quality is " << phase_quality << std::endl;
    }
    
    exit(0);
    
    for (const auto& sp : genotype_posteriors) {
        result.phase_regions.emplace(sp.first, PhaseSet::SamplePhaseRegions {PhaseSet::PhaseRegion {phased_region, 1.0}});
    }
    
    return result;
}

// private methods

    
// non-member methods
    
    bool are_phase_complements(const Genotype<Haplotype>& lhs, const Genotype<Haplotype>& rhs,
                               const std::vector<Variant>& variants)
    {
        if (lhs.ploidy() == 1) return false;
        
        auto regions = get_regions(variants);
        
        return std::all_of(std::cbegin(regions), std::cend(regions),
                           [&lhs, &rhs] (const auto& region) {
                               return splice<Haplotype>(lhs, region) == splice<Haplotype>(rhs, region);
                           });
    }
    
    std::vector<std::vector<Genotype<Haplotype>>>
    partition_phase_complements(const std::vector<Genotype<Haplotype>>& genotypes,
                                const std::vector<Variant>& variants)
    {
        std::vector<std::vector<Genotype<Haplotype>>> result {};
        
        if (genotypes.empty()) return result;
        
        result.reserve(genotypes.size());
        
        result.emplace_back(1, genotypes.front());
        
        std::for_each(std::next(std::cbegin(genotypes)), std::cend(genotypes),
                      [&result, &variants] (const auto& genotype) {
                          auto it = std::find_if(std::begin(result), std::end(result),
                                                 [&genotype, &variants] (const auto& complements) {
                                                     return std::any_of(std::cbegin(complements), std::cend(complements),
                                                                        [&genotype, &variants] (const auto& complement) {
                                                                            return are_phase_complements(genotype, complement, variants);
                                                                        });
                                                 });
                          
                          if (it == std::end(result)) {
                              result.emplace_back(1, genotype);
                          } else {
                              it->emplace_back(genotype);
                          }
                      });
        
        return result;
    }
    
    bool are_phase_complements(const Genotype<Haplotype>& lhs, const Genotype<Haplotype>& rhs)
    {
        if (lhs.ploidy() == 1) return false;
        
        //        std::cout << "comparing " << std::endl;
        //        print_variant_alleles(lhs);
        //        std::cout << std::endl;
        //        print_variant_alleles(rhs);
        //        std::cout << std::endl;
        
        if (lhs == rhs) return true;
        
        // TODO: Generalise... currently diploid only
        const auto& h1 = lhs[0];
        const auto& h2 = lhs[1];
        const auto& h3 = rhs[0];
        const auto& h4 = rhs[1];
        
        auto diff1 = h1.difference(h2);
        
        if (diff1.empty()) return false;
        
        auto diff2 = h3.difference(h4);
        
        if (diff2.empty() || diff1.size() != diff2.size()) return false;
        
        auto regions = get_covered_regions(diff1);
        
        if (regions != get_covered_regions(diff2)) return false;
        
        //        for (const auto& region : regions) {
        //            std::cout << region << std::endl;
        //        }
        
        return std::all_of(std::cbegin(regions), std::cend(regions),
                           [&lhs, &rhs] (const auto& region) {
                               return splice<Haplotype>(lhs, region) == splice<Haplotype>(rhs, region);
                           });
    }
    
    std::vector<std::vector<Genotype<Haplotype>>>
    partition_phase_complements(const std::vector<Genotype<Haplotype>>& genotypes)
    {
        std::vector<std::vector<Genotype<Haplotype>>> result {};
        
        if (genotypes.empty()) return result;
        
        result.reserve(genotypes.size());
        
        result.emplace_back(1, genotypes.front());
        
        std::for_each(std::next(std::cbegin(genotypes)), std::cend(genotypes),
                      [&result] (const auto& genotype) {
                          auto it = std::find_if(std::begin(result), std::end(result),
                                                 [&genotype] (const auto& complements) {
                                                     return std::any_of(std::cbegin(complements), std::cend(complements),
                                                                        [&genotype] (const auto& complement) {
                                                                            return are_phase_complements(genotype, complement);
                                                                        });
                                                 });
                          
                          if (it == std::end(result)) {
                              result.emplace_back(1, genotype);
                          } else {
                              it->emplace_back(genotype);
                          }
                      });
        
        return result;
    }

} // namespace Octopus
