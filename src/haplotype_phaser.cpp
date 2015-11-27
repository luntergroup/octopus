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

bool are_phase_complements(const Genotype<Haplotype>& lhs, const Genotype<Haplotype>& rhs);

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
    GenomeWalker::IndicatorLimit::SharedWithPreviousRegion, GenomeWalker::ExtensionLimit::SharedWithFrontier},
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

HaplotypePhaser::PhaseSets
HaplotypePhaser::phase(const std::vector<Haplotype>& haplotypes,
                       const GenotypePosteriors& genotype_posteriors)
{
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
        
        std::cout << "removing " << low_posterior_haplotypes.size() << " haplotypes" << std::endl;
        
        low_posterior_haplotypes.erase(std::unique(std::begin(low_posterior_haplotypes),
                                                   std::end(low_posterior_haplotypes)),
                                       std::end(low_posterior_haplotypes));
        
        for (const auto& haplotype : low_posterior_haplotypes) {
            tree_.prune_all(haplotype);
        }
    } else {
        tree_.clear();
    }
    
    std::cout << "there are " << buffered_candidates_.size() << " candidates" << std::endl;
    
    std::cout << "current region is " << tree_region_ << std::endl;
    
    next_region_ = walker_.walk(tree_region_, *reads_, buffered_candidates_);
    
    std::cout << "next region will be " << next_region_ << std::endl;
    
    auto phased_region = get_left_overhang(tree_region_, next_region_);
    
    std::cout << "phased region is " << phased_region << std::endl;
    
    tree_.clear(phased_region);
    buffered_candidates_.erase_contained(phased_region);
    
    std::cout << "there are " << buffered_candidates_.size() << " candidates" << std::endl;
    
    exit(0);
    
    PhaseSets result {};
    
    
    
    return result;
}

// private methods

// non-member methods

bool are_phase_complements(const Genotype<Haplotype>& lhs, const Genotype<Haplotype>& rhs)
{
    return true;
}

} // namespace Octopus
