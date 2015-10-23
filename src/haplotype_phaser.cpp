//
//  haplotype_phaser.cpp
//  Octopus
//
//  Created by Daniel Cooke on 30/04/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "haplotype_phaser.hpp"

#include <algorithm>

#include "allele.hpp"
#include "mappable_algorithms.hpp"
#include "search_regions.hpp"

#include <iostream> // DEBUG

//namespace Octopus
//{

// public methods

HaplotypePhaser::HaplotypePhaser(ReferenceGenome& reference)
:
tree_ {reference},
walker_ {max_indicators_, 4, GenomeWalker::IndicatorLimit::NoLimit, GenomeWalker::ExtensionLimit::WithinReadLengthOfFirstIncluded}
{}

HaplotypePhaser::HaplotypePhaser(ReferenceGenome& reference, unsigned max_haplotypes)
:
tree_ {reference},
max_haplotypes_ {max_haplotypes},
walker_ {max_indicators_, 4, GenomeWalker::IndicatorLimit::NoLimit, GenomeWalker::ExtensionLimit::WithinReadLengthOfFirstIncluded}
{}

HaplotypePhaser::HaplotypePhaser(ReferenceGenome& reference, unsigned max_haplotypes, unsigned max_indicators)
:
tree_ {reference},
max_haplotypes_ {max_haplotypes},
max_indicators_ {max_indicators},
walker_ {max_indicators_, 4, GenomeWalker::IndicatorLimit::NoLimit, GenomeWalker::ExtensionLimit::WithinReadLengthOfFirstIncluded}
{}

void HaplotypePhaser::setup(const std::vector<Variant>& candidates, const ReadMap& reads)
{
    if (candidates.empty()) return;
    
    buffered_candidates_.insert(std::cbegin(candidates), std::cend(candidates));
    
    tree_region_ = shift(get_head(buffered_candidates_.leftmost()), -1);
    
    extend_tree(reads);
}

bool HaplotypePhaser::expended_candidates() const noexcept
{
    return buffered_candidates_.empty();
}

std::vector<Haplotype> HaplotypePhaser::get_haplotypes() const
{
    auto haplotypes = tree_.get_haplotypes(tree_region_);
    
    unique_least_complex(haplotypes);
    
    return haplotypes;
}

std::unordered_map<Haplotype, double>
compute_haplotype_posteriors(const std::vector<Haplotype>& haplotypes,
                             const HaplotypePhaser::UnphasedSampleGenotypePosteriors& genotype_posteriors)
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

HaplotypePhaser::PhasedGenotypePosteriors
HaplotypePhaser::phase(const std::vector<Haplotype>& haplotypes,
                       const UnphasedGenotypePosteriors& genotype_posteriors,
                       const ReadMap& reads)
{
    if (max_indicators_ == 0) { // no phasing
        tree_.clear();
        // TODO: need to add an erase method to MappableSet API to rmeove contained range
        auto contained = bases(buffered_candidates_.contained_range(tree_region_));
        buffered_candidates_.erase(contained.begin(), contained.end());
        if (!buffered_candidates_.empty()) {
            tree_region_ = shift(get_head(buffered_candidates_.leftmost()), -1);
        }
    } else {
        std::deque<Haplotype> low_posterior_haplotypes {};
        
        for (const auto& sample_genotype_posteriors : genotype_posteriors) {
            auto sample_haplotype_posteriors = compute_haplotype_posteriors(haplotypes, sample_genotype_posteriors.second);
            
            for (const auto& haplotype_posterior : sample_haplotype_posteriors) {
                if (haplotype_posterior.second < 0.01) {
                    low_posterior_haplotypes.push_back(haplotype_posterior.first);
                }
            }
        }
        
        std::sort(std::begin(low_posterior_haplotypes), std::end(low_posterior_haplotypes));
        
        low_posterior_haplotypes.erase(std::unique(std::begin(low_posterior_haplotypes), std::end(low_posterior_haplotypes)),
                                       std::end(low_posterior_haplotypes));
        
        for (const auto& haplotype : low_posterior_haplotypes) {
            tree_.prune_all(haplotype);
        }
    }
    
    extend_tree(reads);
    
    return genotype_posteriors;
    
//    HaplotypePhaser::PhasedGenotypePosteriors result {};
//    result.reserve(genotype_posteriors.size());
//    
//    return result;
}

// private methods

void extend(HaplotypeTree& tree, const Variant& variant)
{
    tree.extend(variant.get_reference_allele());
    tree.extend(variant.get_alternative_allele());
}

void HaplotypePhaser::extend_tree(const ReadMap& reads)
{
    auto region = walker_.continue_walk(tree_region_, reads, buffered_candidates_);
    
    if (tree_.empty()) {
        tree_region_ = region;
    } else {
        tree_region_ = get_encompassing(tree_region_, region);
    }
    
    auto contained_candidates = buffered_candidates_.contained_range(region);
    
    std::for_each(std::cbegin(contained_candidates), std::cend(contained_candidates),
                  [this] (const auto& candidate) {
                      extend(tree_, candidate);
                  });
}

//} // namespace Octopus
