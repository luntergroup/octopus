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
#include "logging.hpp"

#include <iostream> // DEBUG
#include "timers.hpp"

namespace Octopus
{
auto max_included(const unsigned max_haplotypes)
{
    return 2 * static_cast<unsigned>(std::max(1.0, std::log2(max_haplotypes))) - 1;
}

namespace
{
    auto variants_to_alleles(const MappableFlatSet<Variant>& variants)
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
                                       const MappableFlatSet<Variant>& candidates, const ReadMap& reads,
                                       unsigned max_haplotypes, bool allow_lagging)
:
tree_ {window.get_contig_name(), reference},
walker_ {max_included(max_haplotypes)},
lagged_walker_ {},
alleles_ {variants_to_alleles(candidates)},
reads_ {reads},
current_active_region_ {shift(head_region(alleles_.leftmost(), 0), -1)},
next_active_region_ {},
soft_max_haplotypes_ {max_haplotypes},
rightmost_allele_ {rightmost_mappable(alleles_)},
holdout_set_ {}
{
    if (allow_lagging) {
        lagged_walker_ = GenomeWalker {max_included(max_haplotypes),
                            GenomeWalker::IndicatorLimit::SharedWithPreviousRegion};
    }
    
    if (alleles_.empty()) {
        current_active_region_ = window;
        next_active_region_    = window;
    }
}

const GenomicRegion& HaplotypeGenerator::tell_next_active_region() const
{
    update_next_active_region();
    assert(current_active_region_ <= *next_active_region_);
    return *next_active_region_;
}

std::pair<std::vector<Haplotype>, GenomicRegion> HaplotypeGenerator::progress()
{
    if (alleles_.empty()) {
        return std::make_pair(std::vector<Haplotype> {}, current_active_region_);
    }
    
    holdout_set_.clear(); // TODO: reintroduce holdout alleles and backtrack
    
    update_next_active_region();
    
    assert(next_active_region_ != boost::none);
    
    if (is_after(*next_active_region_, *rightmost_allele_)) {
        return std::make_pair(std::vector<Haplotype> {}, *next_active_region_);
    }
    
    force_forward(*next_active_region_);
    
    auto novel_active_region = *next_active_region_;
    
    if (!tree_.empty()) {
        novel_active_region = right_overhang_region(*next_active_region_, current_active_region_);
    }
    
    current_active_region_ = std::move(*next_active_region_);
    next_active_region_    = boost::none;
    
    const auto novel_active_alleles = overlap_range(alleles_, novel_active_region);
    
    if (novel_active_alleles.empty()) {
        return std::make_pair(std::vector<Haplotype> {}, *next_active_region_);
    }
    
    const auto it = extend_tree_until(novel_active_alleles, tree_, hard_max_haplotypes_);
    
    if (it != std::cend(novel_active_alleles)) {
        assert(holdout_set_.empty());
        
        tree_.remove_overlapped(novel_active_region); // revert
        
        holdout_set_ = compute_holdout_set(current_active_region_);
        
        alleles_.erase_all(std::cbegin(holdout_set_), std::cend(holdout_set_));
        
        Logging::WarningLogger wlog {};
        stream(wlog) << "Holding out " << holdout_set_.size() << " alleles from region "
                     << current_active_region_;
        
        current_active_region_ = walker_.walk(head_region(current_active_region_), reads_, alleles_);
        
        extend_tree(overlap_range(alleles_, current_active_region_), tree_);
    }
    
    auto haplotypes = tree_.extract_haplotypes(calculate_haplotype_region());
    
    if (!is_lagged()) {
        tree_.clear();
    }
    
    return std::make_pair(std::move(haplotypes), current_active_region_);
}

void HaplotypeGenerator::clear_progress() noexcept
{
    tree_.clear();
    holdout_set_.clear();
    next_active_region_ = boost::none;
}

void HaplotypeGenerator::uniquely_keep(const std::vector<Haplotype>& haplotypes)
{
    next_active_region_ = boost::none;
    
    if (!is_active_region_lagged() || tree_.num_haplotypes() == haplotypes.size()) {
        return;
    }
    
    prune_unique(haplotypes, tree_);
}

void HaplotypeGenerator::remove(const std::vector<Haplotype>& haplotypes)
{
    if (haplotypes.empty()) return;
    
    next_active_region_ = boost::none;
    
    if (TRACE_MODE) {
        //Logging::TraceLogger log {};
        //stream(log) << "Removing " << haplotypes.size() << " haplotypes:";
        //debug::print_haplotypes(stream(trace_log), removable_haplotypes);
    } else if (DEBUG_MODE) {
        Logging::DebugLogger log {};
        stream(log) << "Discarding " << haplotypes.size() << " haplotypes from generator";
    }
    
    if (!is_active_region_lagged() || haplotypes.size() == tree_.num_haplotypes()) {
        tree_.clear();
        alleles_.erase_overlapped(current_active_region_);
    } else {
        prune_all(haplotypes, tree_);
    }
}

void HaplotypeGenerator::remove(const std::vector<std::reference_wrapper<const Haplotype>>& haplotypes)
{
    if (haplotypes.empty()) return;
    
    next_active_region_ = boost::none;
    
    if (TRACE_MODE) {
        //Logging::TraceLogger log {};
        //stream(log) << "Removing " << haplotypes.size() << " haplotypes:";
        //debug::print_haplotypes(stream(trace_log), removable_haplotypes);
    } else if (DEBUG_MODE) {
        Logging::DebugLogger log {};
        stream(log) << "Discarding " << haplotypes.size() << " haplotypes from generator";
    }
    
    if (!is_active_region_lagged() || haplotypes.size() == tree_.num_haplotypes()) {
        tree_.clear();
        alleles_.erase_overlapped(current_active_region_);
    } else {
        prune_all(haplotypes, tree_);
    }
}

template <typename Range>
bool can_remove_entire_passed_region(const GenomicRegion& current_active_region,
                                     const GenomicRegion& next_active_region,
                                     const Range& passed_alleles)
{
    return passed_alleles.empty()
            || (is_after(next_active_region, current_active_region)
                && !overlaps(rightmost_region(passed_alleles), next_active_region));
}

template <typename Range>
bool requires_staged_removal(const Range& passed_alleles)
{
    if (passed_alleles.empty()) return false;
    
    const auto it = std::find_if_not(std::next(std::crbegin(passed_alleles)),
                                     std::crend(passed_alleles),
                                     [&passed_alleles] (const auto& allele) {
                                         return overlaps(allele, passed_alleles.back());
                                     });
    
    return it == std::crend(passed_alleles) || is_position(*it);
}

void HaplotypeGenerator::force_forward(GenomicRegion to)
{
    assert(to > current_active_region_);
    
    next_active_region_ = std::move(to);
    
    if (begins_before(current_active_region_, *next_active_region_)) {
        auto passed_region = left_overhang_region(current_active_region_, *next_active_region_);
        
        const auto passed_alleles = overlap_range(alleles_, passed_region);
        
        if (passed_alleles.empty()) return;
        
        if (can_remove_entire_passed_region(current_active_region_, *next_active_region_,
                                            passed_alleles)) {
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
    } else if (is_after(*next_active_region_, current_active_region_)) {
        tree_.clear();
    }
}

// private methods

bool HaplotypeGenerator::is_lagged() const noexcept
{
    return lagged_walker_ != boost::none;
}

bool HaplotypeGenerator::is_active_region_lagged() const
{
    if (!is_lagged()) return false;
    
    const auto next_lagged_region = lagged_walker_->walk(current_active_region_, reads_, alleles_);
    
    return overlaps(current_active_region_, next_lagged_region);
}

void HaplotypeGenerator::update_next_active_region() const
{
    if (!next_active_region_) {
        if (is_lagged()) {
            if (contains(current_active_region_, *rightmost_allele_)) {
                next_active_region_ = shift(tail_region(*rightmost_allele_), 2);
                return;
            }
            
            auto max_lagged_region = lagged_walker_->walk(current_active_region_, reads_,
                                                          alleles_);
            
            if (!overlaps(current_active_region_, max_lagged_region)) {
                next_active_region_ = std::move(max_lagged_region);
            } else {
                auto tmp_tree = tree_;
                
                if (begins_before(current_active_region_, max_lagged_region)) {
                    const auto passed_region  = left_overhang_region(current_active_region_,
                                                                     max_lagged_region);
                    const auto passed_alleles = overlap_range(alleles_, passed_region);
                    
                    if (requires_staged_removal(passed_alleles)) {
                        const auto first_removal_region = expand_rhs(passed_region, -1);
                        tmp_tree.remove_overlapped(first_removal_region);
                        tmp_tree.remove_overlapped(tail_region(first_removal_region));
                    } else {
                        tmp_tree.remove_overlapped(expand_rhs(passed_region, -1));
                    }
                }
                
                const auto indicator_region = overlapped_region(current_active_region_,
                                                                max_lagged_region);
                const auto indicator_alleles = overlap_range(alleles_, indicator_region);
                
                assert(!indicator_alleles.empty());
                
                const auto mutually_exclusive_indicator_regions = extract_covered_regions(indicator_alleles);
                
                for (const auto& region : mutually_exclusive_indicator_regions) {
                    if (tmp_tree.num_haplotypes() < soft_max_haplotypes_) {
                        break;
                    }
                    tmp_tree.remove_overlapped(region);
                }
                
                const auto novel_region  = right_overhang_region(max_lagged_region,
                                                               current_active_region_);
                const auto novel_alleles = overlap_range(alleles_, novel_region);
                
                assert(!novel_alleles.empty());
                
                const auto mutually_exclusive_novel_regions = extract_covered_regions(novel_alleles);
                
                unsigned num_regions_added {0};
                for (const auto& region : mutually_exclusive_novel_regions) {
                    const auto interacting_alleles = overlap_range(novel_alleles, region);
                    
                    assert(!interacting_alleles.empty());
                    
                    const auto it = extend_tree_until(interacting_alleles, tmp_tree,
                                                      hard_max_haplotypes_);
                    
                    if (it != std::cend(interacting_alleles)) {
                        next_active_region_ = novel_region; // Revert to non-lagged behaviour
                        return;
                    }
                    
                    ++num_regions_added;
                    
                    if (tmp_tree.num_haplotypes() >= soft_max_haplotypes_) {
                        break;
                    }
                }
                
                next_active_region_ = tmp_tree.encompassing_region();
                
                if (num_regions_added > 0
                    && is_empty(mutually_exclusive_novel_regions[num_regions_added - 1])) {
                    // To ensure we progress on a single novel insertion
                    next_active_region_ = expand_rhs(*next_active_region_, 1);
                }
                
                if (*next_active_region_ == current_active_region_) {
                    assert(num_regions_added == 0);
                    next_active_region_ = walker_.walk(current_active_region_, reads_, alleles_);
                } else if (begins_before(current_active_region_, *next_active_region_)
                           && is_empty_region(indicator_alleles.front())) {
                    // To be clearer about insertion containment
                    next_active_region_ = expand_lhs(*next_active_region_, 1);
                }
            }
        } else {
            next_active_region_ = walker_.walk(current_active_region_, reads_, alleles_);
        }
    }
}

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
        
        if (std::exp2(count_overlapped(cur, re)) < hard_max_haplotypes_) break;
    }
    
    return MappableFlatMultiSet<Allele> {it, std::end(tmp)};
}

template <typename Range>
auto sum_indel_sizes(const Range& alleles)
{
    return std::accumulate(std::cbegin(alleles), std::cend(alleles),
                           GenomicRegion::SizeType {0},
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
    const auto overlapped = overlap_range(alleles_, current_active_region_);
    
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
