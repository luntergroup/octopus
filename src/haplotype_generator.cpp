//
//  haplotype_generator.cpp
//  Octopus
//
//  Created by Daniel Cooke on 20/02/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#include "haplotype_generator.hpp"

#include <algorithm>
#include <deque>
#include <iterator>
#include <cassert>

#include "variant.hpp"

#include <iostream> // DEBUG

namespace Octopus
{
    auto max_indcluded(const unsigned max_haplotypes)
    {
        return 2 * static_cast<unsigned>(std::max(1.0, std::log2(max_haplotypes)));
    }
    
    namespace
    {
        auto variants_to_alleles(const std::vector<Variant>& variants)
        {
            std::deque<Allele> alleles {};
            
            for (const auto& variant : variants) {
                alleles.push_back(variant.get_ref_allele());
                alleles.push_back(variant.get_alt_allele());
            }
            
            std::sort(std::begin(alleles), std::end(alleles));
            
            const auto it = std::unique(std::begin(alleles), std::end(alleles));
            
            const auto it1 = std::make_move_iterator(std::begin(alleles));
            const auto it2 = std::make_move_iterator(it);
            
            return MappableSet<Allele> {it1, it2};
        }
    }
    
    // public members
    
    HaplotypeGenerator::HaplotypeGenerator(const GenomicRegion& window, const ReferenceGenome& reference,
                                           const std::vector<Variant>& candidates, const ReadMap& reads,
                                           unsigned max_haplotypes, unsigned max_indicators)
    :
    tree_ {window.get_contig_name(), reference},
    walker_ {2 * max_indicators, max_indcluded(max_haplotypes)},
    alleles_ {variants_to_alleles(candidates)},
    reads_ {reads},
    current_active_region_ {shift(head_region(alleles_.leftmost(), 0), -1)},
    next_active_region_ {},
    max_haplotypes_ {max_haplotypes},
    active_allele_counts_ {}
    {}
    
    GenomicRegion HaplotypeGenerator::tell_next_active_region()
    {
        if (!next_active_region_) {
            next_active_region_ = walker_.walk(current_active_region_, reads_, alleles_);
        }
        
        assert(current_active_region_ < *next_active_region_);
        
        return *next_active_region_;
    }
    
    namespace
    {
        auto count_contained(const Allele& allele, const std::vector<Haplotype>& haplotypes)
        {
            return std::count_if(std::cbegin(haplotypes), std::cend(haplotypes),
                                 [&allele] (const auto& haplotype) {
                                     return haplotype.contains_exact(allele);
                                 });
        }
    } // namespace
    
    std::pair<std::vector<Haplotype>, GenomicRegion> HaplotypeGenerator::progress()
    {
        if (!next_active_region_) {
            next_active_region_ = walker_.walk(current_active_region_, reads_, alleles_);
        }
        
        force_forward(*next_active_region_);
        
        extend_tree(alleles_.overlap_range(*next_active_region_), tree_);
        
        current_active_region_ = std::move(*next_active_region_);
        
        next_active_region_ = boost::none;
        
        auto haplotype_region = calculate_haplotype_region();
        
        assert(contains(haplotype_region, current_active_region_));
        
        auto haplotypes = tree_.extract_haplotypes(haplotype_region);
        
        active_allele_counts_.clear();
        
        for (const auto& allele : alleles_.overlap_range(current_active_region_)) {
            active_allele_counts_.emplace(allele, count_contained(allele, haplotypes));
        }
        
        return std::make_pair(std::move(haplotypes), current_active_region_);
    }
    
    bool is_not_in_haplotypes(const Allele& allele, const std::vector<Haplotype>& haplotypes)
    {
        return std::none_of(std::cbegin(haplotypes), std::cend(haplotypes),
                            [&allele] (const auto& haplotype) {
                                return haplotype.contains(allele);
                            });
    }
    
    void HaplotypeGenerator::keep_haplotypes(const std::vector<Haplotype>& haplotypes)
    {
        for (const auto& haplotype : haplotypes) {
            tree_.prune_unique(haplotype);
        }
        
        std::deque<Allele> unused_alleles {};
        
        for (const auto& candidate : alleles_.overlap_range(current_active_region_)) {
            if (is_not_in_haplotypes(candidate, haplotypes)) {
                unused_alleles.push_back(candidate);
            }
        }
        
        for (const auto& allele : unused_alleles) {
            alleles_.erase(allele);
        }
        
        active_allele_counts_.clear();
        
        for (const auto& allele : alleles_.overlap_range(current_active_region_)) {
            active_allele_counts_.emplace(allele, count_contained(allele, haplotypes));
        }
    }
    
    void HaplotypeGenerator::remove_haplotypes(const std::vector<Haplotype>& haplotypes)
    {
        for (const auto& haplotype : haplotypes) {
            tree_.prune_all(haplotype);
        }
        
        std::deque<Allele> removed_alleles {};
        
        for (const auto& allele : alleles_.overlap_range(current_active_region_)) {
            if (active_allele_counts_[allele] == count_contained(allele, haplotypes)) {
                removed_alleles.push_back(allele);
            }
        }
        
        for (const auto& allele : removed_alleles) {
            alleles_.erase(allele);
            active_allele_counts_.erase(allele);
        }
    }
    
    void HaplotypeGenerator::force_forward(GenomicRegion to)
    {
        assert(to > current_active_region_);
        
        next_active_region_ = std::move(to);
        
        if (begins_before(current_active_region_, *next_active_region_)) {
            const auto passed_region = left_overhang_region(current_active_region_, *next_active_region_);
            alleles_.erase_overlapped(passed_region);
            tree_.remove(passed_region);
        }
    }
    
    // private methods
    
    GenomicRegion HaplotypeGenerator::calculate_haplotype_region() const
    {
        // TODO
        return expand(current_active_region_, 100);
    }
    
} // namespace Octopus
