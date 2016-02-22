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
        return static_cast<unsigned>(std::max(1.0, std::log2(max_haplotypes)));
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
    walker_ {max_indicators, max_indcluded(max_haplotypes)},
    alleles_ {variants_to_alleles(candidates)},
    reads_ {reads},
    current_region_ {shift(head_region(alleles_.leftmost(), 0), -1)},
    next_region_ {},
    max_haplotypes_ {max_haplotypes}
    {}
    
    bool HaplotypeGenerator::done() const noexcept
    {
        return alleles_.empty();
    }
    
    GenomicRegion HaplotypeGenerator::tell_next_active_region()
    {
        if (!next_region_) {
            next_region_ = walker_.walk(current_region_, reads_, alleles_);
        }
        
        assert(ends_before(current_region_, *next_region_));
        
        return *next_region_;
    }
    
    std::pair<std::vector<Haplotype>, GenomicRegion> HaplotypeGenerator::progress()
    {
        if (!next_region_) {
            next_region_ = walker_.walk(current_region_, reads_, alleles_);
        }
        
        assert(ends_before(current_region_, *next_region_));
        
        if (begins_before(current_region_, *next_region_)) {
            const auto passed_region = left_overhang_region(current_region_, *next_region_);
            alleles_.erase_contained(passed_region);
            tree_.remove(passed_region);
        }
        
        extend_tree(alleles_.overlap_range(*next_region_), tree_);
        
        const auto haplotype_region = expand(current_region_, 50);
        
        current_region_ = std::move(*next_region_);
        
        next_region_ = boost::none;
        
        return std::make_pair(tree_.extract_haplotypes(haplotype_region), current_region_);
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
        
        std::deque<std::reference_wrapper<const Allele>> unused_alleles {};
        
        for (const auto& candidate : alleles_.overlap_range(current_region_)) {
            if (is_not_in_haplotypes(candidate, haplotypes)) {
                unused_alleles.push_back(candidate);
            }
        }
        
        for (const auto& allele : unused_alleles) {
            alleles_.erase(allele);
        }
    }
    
    void HaplotypeGenerator::force_forward(GenomicRegion to)
    {
        assert(to > current_region_);
        next_region_ = std::move(to);
    }
    
} // namespace Octopus
