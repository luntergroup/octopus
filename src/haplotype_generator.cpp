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
#include <numeric>

#include "variant.hpp"
#include "haplotype.hpp"
#include "mappable.hpp"
#include "mappable_algorithms.hpp"

#include <iostream> // DEBUG
#include "timers.hpp"

namespace Octopus
{
    auto max_indcluded(const unsigned max_haplotypes)
    {
        return 2 * static_cast<unsigned>(std::max(1.0, std::log2(max_haplotypes))) - 1;
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
            
            return MappableFlatSet<Allele> {it1, it2};
        }
    }
    
    namespace debug
    {
        template <typename Range>
        void print_active_alleles(const Range& alleles, const GenomicRegion& active_region);
    } // namespace debug
    
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
    is_lagged_ {max_indicators > 0},
    holdout_set_ {},
    active_allele_counts_ {}
    {}
    
    GenomicRegion HaplotypeGenerator::tell_next_active_region()
    {
        if (!next_active_region_) {
            next_active_region_ = walker_.walk(current_active_region_, reads_, alleles_);
        }
        
        assert(current_active_region_ <= *next_active_region_);
        
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
    
    template <typename Range>
    auto estimate_num_haplotypes(const Range& alleles)
    {
        return std::exp2(size(alleles) / 2);
    }
    
    std::pair<std::vector<Haplotype>, GenomicRegion> HaplotypeGenerator::progress()
    {
        holdout_set_.clear(); // TODO: reintroduce holdout alleles and backtrack
        
        if (!next_active_region_) {
            next_active_region_ = walker_.walk(current_active_region_, reads_, alleles_);
        }
        
        if (*next_active_region_ == current_active_region_) {
            return std::make_pair(std::vector<Haplotype> {}, current_active_region_);
        }
        
        force_forward(*next_active_region_);
        
        current_active_region_ = std::move(*next_active_region_);
        next_active_region_    = boost::none;
        
        auto active_alleles = alleles_.overlap_range(current_active_region_);
        
        if (estimate_num_haplotypes(active_alleles) > hard_max_haplotypes_) {
            assert(holdout_set_.empty());
            
            holdout_set_ = compute_holdout_set(current_active_region_);
            
            alleles_.erase_all(std::cbegin(holdout_set_), std::cend(holdout_set_));
            
            std::cout << "Warning: holding out " << holdout_set_.size() << " alleles from region "
                        << current_active_region_ << '\n';
            
            current_active_region_ = walker_.walk(head_region(current_active_region_), reads_, alleles_);
            active_alleles = alleles_.overlap_range(current_active_region_);
        }
        
        //debug::print_active_alleles(active_alleles, current_active_region_);
        
        extend_tree(active_alleles, tree_);
        
        auto haplotypes = tree_.extract_haplotypes(calculate_haplotype_region());
        
        if (!is_lagged_) {
            tree_.clear();
        }
        
        active_allele_counts_.clear();
        
        for (const auto& allele : active_alleles) {
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
        if (!is_lagged_ || tree_.num_haplotypes() == haplotypes.size()) return;
        
        for (const auto& haplotype : haplotypes) {
            tree_.prune_unique(haplotype);
        }
        
        std::deque<Allele> unused_alleles {};
        
        for (const auto& candidate : alleles_.overlap_range(current_active_region_)) {
            if (is_not_in_haplotypes(candidate, haplotypes)) {
                unused_alleles.push_back(candidate);
            }
        }
        
        alleles_.erase_all(std::cbegin(unused_alleles), std::cend(unused_alleles));
        
        active_allele_counts_.clear();
        
        for (const auto& allele : alleles_.overlap_range(current_active_region_)) {
            active_allele_counts_.emplace(allele, count_contained(allele, haplotypes));
        }
    }
    
    void HaplotypeGenerator::remove_haplotypes(const std::vector<Haplotype>& haplotypes)
    {
        if (!is_lagged_ || haplotypes.size() == tree_.num_haplotypes()) {
            tree_.clear();
            alleles_.erase_overlapped(current_active_region_);
            active_allele_counts_.clear();
        } else {
            for (const auto& haplotype : haplotypes) {
                tree_.prune_all(haplotype);
            }
            
            std::deque<Allele> removed_alleles {};
            
            for (const auto& allele : alleles_.overlap_range(current_active_region_)) {
                if (active_allele_counts_[allele] == count_contained(allele, haplotypes)) {
                    removed_alleles.push_back(allele);
                }
            }
            
            alleles_.erase_all(std::cbegin(removed_alleles), std::cend(removed_alleles));
            
            for (const auto& allele : removed_alleles) {
                active_allele_counts_.erase(allele);
            }
        }
    }
    
    template <typename Range>
    bool requires_staged_removal(const Range& passed_alleles)
    {
        const auto rbegin = std::make_reverse_iterator(std::prev(std::cend(passed_alleles)));
        const auto rend   = std::make_reverse_iterator(std::cbegin(passed_alleles));
        
        const auto it = std::find_if_not(rbegin, rend,
                                         [&passed_alleles] (const auto& allele) {
                                             return overlaps(allele, passed_alleles.back());
                                         });
        
        return it == rend || is_position(*it);
    }
    
    void HaplotypeGenerator::force_forward(GenomicRegion to)
    {
        assert(to > current_active_region_);
        
        next_active_region_ = std::move(to);
        
        if (begins_before(current_active_region_, *next_active_region_)) {
            auto passed_region = left_overhang_region(current_active_region_, *next_active_region_);
            
            const auto passed_alleles = alleles_.overlap_range(passed_region);
            
            if (passed_alleles.empty()) return;
            
            if (!is_empty_region(passed_alleles.back())) {
                alleles_.erase_overlapped(passed_region);
                tree_.remove_overlapped(passed_region);
            } else if (requires_staged_removal(passed_alleles)) {
                // We need to be careful here as insertions adjacent to passed_region are
                // considered overlapped and would be wrongly erased if we erased the whole
                // region. But, we also want to clear all single base alleles left adjacent with
                // next_active_region_, as they have truly been passed.
                
                // This will erase everthing to the left of the adjacent insertion, other than
                // the single base alleles adjacent with next_active_region_.
                const auto first_removal_region = expand_rhs(passed_region, -1);
                
                alleles_.erase_overlapped(first_removal_region);
                tree_.remove_overlapped(first_removal_region);
                
                // This will erase the remaining single base alleles in passed_region, but not the
                // insertions in next_active_region_.
                const auto second_removal_region = tail_region(first_removal_region);
                
                alleles_.erase_overlapped(second_removal_region);
                tree_.remove_overlapped(second_removal_region);
            } else {
                const auto removal_region = expand_rhs(passed_region, -1);
                
                alleles_.erase_overlapped(removal_region);
                tree_.remove_overlapped(removal_region);
            }
        }
    }
    
    // private methods
    
    MappableFlatMultiSet<Allele>
    HaplotypeGenerator::compute_holdout_set(const GenomicRegion& active_region) const
    {
        auto active_alleles = copy_overlapped(alleles_, active_region);
        
        std::vector<Allele> tmp {std::cbegin(active_alleles), std::cend(active_alleles)};
        
        auto it = std::end(tmp);
        
        while (true) {
            const auto r = largest_region(std::begin(tmp), it);
            
            it = std::partition(std::begin(tmp), it,
                                [&r] (const Allele& allele) {
                                    return !is_same_region(allele, r);
                                });
            
            MappableFlatSet<Allele> cur {std::begin(tmp), it};
            
            auto re = walker_.walk(head_region(cur.leftmost()), reads_, cur);
            
            if (std::exp2(cur.count_overlapped(re)) < hard_max_haplotypes_) break;
        }
        
        return MappableFlatMultiSet<Allele> {it, std::end(tmp)};
    }
    
    template <typename Range>
    auto sum_indel_sizes(const Range& alleles)
    {
        return std::accumulate(std::cbegin(alleles), std::cend(alleles), GenomicRegion::SizeType {0},
                               [] (const auto curr, const Allele& allele) {
                                   if (is_insertion(allele)) {
                                       return curr + sequence_size(allele);
                                   } else if (is_deletion(allele)) {
                                       return curr + region_size(allele);
                                   }
                                   return curr;
                               });
    }
    
    GenomicRegion HaplotypeGenerator::calculate_haplotype_region() const
    {
        const auto overlapped = alleles_.overlap_range(current_active_region_);
        
        // We want to keep haplotypes as small as possible, while allowing sufficient flanking
        // reference sequence for full read re-mapping and alignment (i.e. the read must be
        // contained by the haplotype). Note the sum of the indel sizes may not be sufficient
        // as the candidate generator may not propopse all variation in the original reads.
        const auto additional_padding = 2 * sum_indel_sizes(overlapped) + 30;
        
        const auto& lhs_read = *leftmost_overlapped(reads_.get(), current_active_region_);
        const auto& rhs_read = *rightmost_overlapped(reads_.get(), current_active_region_);
        
        const auto unpadded_region = encompassing_region(lhs_read, rhs_read);
        
        if (region_begin(lhs_read) < additional_padding / 2) {
            const auto lhs_padding = region_begin(lhs_read);
            const auto rhs_padding = additional_padding - lhs_padding;
            return expand_lhs(expand_rhs(unpadded_region, rhs_padding), lhs_padding);
        }
        
        return expand(unpadded_region, additional_padding / 2);
    }
    
    namespace debug
    {
        template <typename Range>
        void print_active_alleles(const Range& alleles, const GenomicRegion& active_region)
        {
            std::cout << "printing " << std::distance(std::cbegin(alleles), std::cend(alleles))
                        << " alleles in active region " << active_region << '\n';
            std::copy(std::cbegin(alleles), std::cend(alleles),
                      std::ostream_iterator<Allele>(std::cout, "\n"));
        }
    } // namespace debug
} // namespace Octopus
