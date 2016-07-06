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
#include <numeric>
#include <cassert>

#include "variant.hpp"
#include "haplotype.hpp"
#include "mappable.hpp"
#include "mappable_algorithms.hpp"
#include "logging.hpp"

#include <iostream> // DEBUG
#include "timers.hpp"

namespace Octopus
{
// HaplotypeOverflowError

HaplotypeGenerator::HaplotypeOverflowError::HaplotypeOverflowError(GenomicRegion region)
:
runtime_error {"HaplotypeOverflowError"},
region_ {std::move(region)},
message_ {}
{}

const char* HaplotypeGenerator::HaplotypeOverflowError::what() const noexcept
{
    return runtime_error::what();
}

const GenomicRegion& HaplotypeGenerator::HaplotypeOverflowError::region() const noexcept
{
    return region_;
}

unsigned HaplotypeGenerator::HaplotypeOverflowError::overflow_size() const noexcept
{
    return 0;
}

// HaplotypeGenerator

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
            alleles.push_back(variant.ref_allele());
            alleles.push_back(variant.alt_allele());
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
                                       HaplotypeLimits haplotype_limits,
                                       LaggingPolicy lagging, Haplotype::SizeType min_pad)
:
tree_ {window.contig_name(), reference},
walker_ {max_included(haplotype_limits.soft_max)},
lagged_walker_ {},
haplotype_limits_ {haplotype_limits},
min_pad_ {min_pad},
alleles_ {variants_to_alleles(candidates)},
reads_ {reads},
current_active_region_ {shift(head_region(alleles_.leftmost(), 0), -1)},
next_active_region_ {},
holdout_set_ {},
current_holdout_region_ {},
previous_holdout_regions_ {}
{
    if (lagging != LaggingPolicy::None) {
        GenomeWalker::IndicatorPolicy policy;
        
        if (lagging == LaggingPolicy::Conservative) {
            policy = GenomeWalker::IndicatorPolicy::IncludeIfSharedWithNovelRegion;
        } else {
            policy = GenomeWalker::IndicatorPolicy::IncludeIfLinkableToNovelRegion;
        }
        
        lagged_walker_ = GenomeWalker {max_included(haplotype_limits_.soft_max), policy};
    }
    
    if (alleles_.empty()) {
        current_active_region_ = window;
        next_active_region_    = window;
    } else {
        rightmost_allele_ = *rightmost_mappable(alleles_);
    }
}

GenomicRegion HaplotypeGenerator::tell_next_active_region() const
{
    update_next_active_region();
    assert(current_active_region_ <= *next_active_region_ || !holdout_set_.empty());
    return *next_active_region_;
}

template <typename Range>
bool can_remove_entire_passed_region(const GenomicRegion& current_active_region,
                                     const GenomicRegion& next_active_region,
                                     const Range& passed_alleles)
{
    return passed_alleles.empty() || !overlaps(rightmost_region(passed_alleles), next_active_region);
}

template <typename Range>
bool requires_staged_removal(const Range& passed_alleles)
{
    if (passed_alleles.empty() || !is_empty_region(passed_alleles.back())) {
        return false;
    }
    
    auto it = std::find_if_not(std::next(std::crbegin(passed_alleles)), std::crend(passed_alleles),
                               [&passed_alleles] (const auto& allele) {
                                   return overlaps(allele, passed_alleles.back());
                               });
    
    return it == std::crend(passed_alleles) || is_position(*std::prev(it));
}

void HaplotypeGenerator::progress(GenomicRegion to)
{
    if (to == current_active_region_) return;
    
    if (to < current_active_region_ && holdout_set_.empty()) {
        throw std::runtime_error {"HaplotypeGenerator: trying to progress backwards from "
            + to_string(current_active_region_) + " to " + to_string(to)};
    }
    
    next_active_region_ = std::move(to);
    
    if (holdout_set_.empty() && begins_before(current_active_region_, *next_active_region_)) {
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

void HaplotypeGenerator::stop() noexcept
{
    tree_.clear();
    reset_next_active_region();
}

HaplotypeGenerator::HaplotypePacket HaplotypeGenerator::generate()
{
    if (alleles_.empty()) {
        return std::make_tuple(std::vector<Haplotype> {}, current_active_region_, true);
    }
    
    if (try_reintroducing_holdout_set()) {
        next_active_region_ = boost::none;
        
        if (tree_.num_haplotypes() > haplotype_limits_.hard_max) {
            throw HaplotypeOverflowError {current_active_region_};
        }
    } else {
        update_next_active_region();
        
        assert(next_active_region_ != boost::none);
        
        if (is_after(*next_active_region_, *rightmost_allele_)) {
            return std::make_tuple(std::vector<Haplotype> {}, *next_active_region_, true);
        }
        
        progress(*next_active_region_);
        
        auto novel_active_region = *next_active_region_;
        
        if (!tree_.empty()) {
            novel_active_region = right_overhang_region(*next_active_region_, current_active_region_);
        }
        
        current_active_region_ = *std::move(next_active_region_);
        
        reset_next_active_region();
        
        auto novel_active_alleles = overlap_range(alleles_, novel_active_region);
        
        if (novel_active_alleles.empty()) {
            return std::make_tuple(std::vector<Haplotype> {}, current_active_region_, true);
        }
        
        auto it = extend_tree_until(novel_active_alleles, tree_, haplotype_limits_.holdout_max);
        
        if (it != std::cend(novel_active_alleles)) {
            assert(holdout_set_.empty());
            
            tree_.remove_overlapped(novel_active_region); // reverts tree
            
            try {
                set_holdout_set(novel_active_region);
                current_active_region_ = head_region(current_active_region_);
                update_next_active_region();
                progress(*next_active_region_);
                current_active_region_ = *std::move(next_active_region_);
                reset_next_active_region();
                extend_tree_until(overlap_range(alleles_, current_active_region_), tree_,
                                  haplotype_limits_.soft_max);
            } catch (const HaplotypeOverflowError&) {
                it = extend_tree_until(novel_active_alleles, tree_, haplotype_limits_.hard_max);
                
                if (it != std::cend(novel_active_alleles)) {
                    throw;
                }
                
                current_active_region_ = tree_.encompassing_region();
            }
        }
    }
    
    auto haplotypes = tree_.extract_haplotypes(calculate_haplotype_region());
    
    if (!is_lagging_enabled()) tree_.clear();
    
    return std::make_tuple(std::move(haplotypes), current_active_region_, holdout_set_.empty());
}

bool HaplotypeGenerator::removal_has_impact() const
{
    if (!holdout_set_.empty()) return true;
    
    if (!is_lagging_enabled() || contains(current_active_region_, *rightmost_allele_)) return false;
    
    const auto max_lagged_region = lagged_walker_->walk(current_active_region_, reads_, alleles_);
    
    return overlaps(max_lagged_region, current_active_region_);
}

unsigned HaplotypeGenerator::max_removal_impact() const
{
    if (!holdout_set_.empty()) return tree_.num_haplotypes();
    
    if (!is_lagging_enabled() || contains(current_active_region_, *rightmost_allele_)) return 0;
    
    const auto max_lagged_region = lagged_walker_->walk(current_active_region_, reads_, alleles_);
    
    if (!overlaps(max_lagged_region, current_active_region_)) return 0;
    
    const auto novel_region = right_overhang_region(max_lagged_region, current_active_region_);
    
    const auto num_novel_alleles = count_overlapped(alleles_, novel_region);
    
    if (num_novel_alleles == 0) return 0;
    
    const auto max_new_haplotypes = std::max(static_cast<unsigned>(std::exp2(num_novel_alleles / 2)), 1u);
    
    const auto num_leftover_haplotypes = haplotype_limits_.soft_max / max_new_haplotypes;
    
    const auto cur_num_haplotypes = tree_.num_haplotypes();
    
    if (cur_num_haplotypes > num_leftover_haplotypes) {
        return cur_num_haplotypes - num_leftover_haplotypes;
    }
    
    return cur_num_haplotypes;
}

// private methods

bool HaplotypeGenerator::is_lagging_enabled() const noexcept
{
    return lagged_walker_ != boost::none;
}

bool HaplotypeGenerator::is_active_region_lagged() const
{
    if (!is_lagging_enabled()) return false;
    
    const auto next_lagged_region = lagged_walker_->walk(current_active_region_, reads_, alleles_);
    
    return overlaps(current_active_region_, next_lagged_region);
}

void HaplotypeGenerator::reset_next_active_region() const noexcept
{
    next_active_region_ = boost::none;
}

void HaplotypeGenerator::update_next_active_region() const
{
    if (!next_active_region_) {
        if (is_lagging_enabled()) {
            if (contains(current_active_region_, *rightmost_allele_)) {
                next_active_region_ = shift(tail_region(*rightmost_allele_), 2);
                return;
            }
            
            auto max_lagged_region = lagged_walker_->walk(current_active_region_, reads_, alleles_);
            
            if (!overlaps(current_active_region_, max_lagged_region)) {
                next_active_region_ = std::move(max_lagged_region);
            } else {
                auto test_tree = tree_; // use a temporary tree to see how much we can lag
                
                if (begins_before(current_active_region_, max_lagged_region)) {
                    const auto novel_region  = right_overhang_region(max_lagged_region, current_active_region_);
                    const auto novel_alleles = overlap_range(alleles_, novel_region);
                    
                    const auto it = extend_tree_until(novel_alleles, test_tree, haplotype_limits_.soft_max);
                    
                    test_tree.remove_overlapped(novel_region); // undo previous extension
                    
                    if (it == std::cend(novel_alleles)) {
                        max_lagged_region = encompassing_region(current_active_region_, max_lagged_region);
                    } else {
                        const auto passed_region  = left_overhang_region(current_active_region_,
                                                                         max_lagged_region);
                        const auto passed_alleles = overlap_range(alleles_, passed_region);
                        
                        if (can_remove_entire_passed_region(current_active_region_,
                                                            max_lagged_region, passed_alleles)) {
                            test_tree.remove_overlapped(passed_region);
                        } else if (requires_staged_removal(passed_alleles)) {
                            const auto first_removal_region = expand_rhs(passed_region, -1);
                            test_tree.remove_overlapped(first_removal_region);
                            test_tree.remove_overlapped(tail_region(first_removal_region));
                        } else {
                            test_tree.remove_overlapped(expand_rhs(passed_region, -1));
                        }
                    }
                }
                
                const auto indicator_region = overlapped_region(current_active_region_, max_lagged_region);
                const auto indicator_alleles = overlap_range(alleles_, indicator_region);
                
                assert(!indicator_alleles.empty());
                
                const auto mutually_exclusive_indicator_regions = extract_mutually_exclusive_regions(indicator_alleles);
                
                for (const auto& region : mutually_exclusive_indicator_regions) {
                    if (test_tree.num_haplotypes() < haplotype_limits_.soft_max) {
                        break;
                    }
                    test_tree.remove_overlapped(region);
                }
                
                const auto novel_region  = right_overhang_region(max_lagged_region, current_active_region_);
                const auto novel_alleles = overlap_range(alleles_, novel_region);
                
                assert(!novel_alleles.empty());
                
                auto mutually_exclusive_novel_regions = extract_mutually_exclusive_regions(novel_alleles);
                
                if (mutually_exclusive_indicator_regions.back() == mutually_exclusive_novel_regions.front()) {
                    assert(is_empty(mutually_exclusive_novel_regions.front()));
                    // always consider lhs boundry insertions
                    mutually_exclusive_novel_regions.erase(std::begin(mutually_exclusive_novel_regions));
                }
                
                assert(!mutually_exclusive_indicator_regions.empty());
                
                for (const auto& region : mutually_exclusive_novel_regions) {
                    const auto interacting_alleles = overlap_range(novel_alleles, region);
                    
                    const auto it = extend_tree_until(interacting_alleles, test_tree,
                                                      haplotype_limits_.hard_max);
                    
                    if (it != std::cend(interacting_alleles)) {
                        test_tree.clear();
                        break;
                    }
                    
                    if (test_tree.num_haplotypes() >= haplotype_limits_.soft_max) {
                        if (region != mutually_exclusive_novel_regions.front()) {
                            test_tree.remove_overlapped(region); // always add one novel region
                        }
                        break;
                    }
                }
                
                if (!test_tree.empty()) {
                    next_active_region_ = test_tree.encompassing_region();
                } else {
                    next_active_region_ = novel_region; // revert to non-lagged behaviour
                }
                
                if (*next_active_region_ == current_active_region_) {
                    next_active_region_ = walker_.walk(current_active_region_, reads_, alleles_);
                }
            }
        } else {
            next_active_region_ = walker_.walk(current_active_region_, reads_, alleles_);
        }
    }
}

void HaplotypeGenerator::set_holdout_set(const GenomicRegion& active_region)
{
    assert(holdout_set_.empty());
    
    if (previous_holdout_regions_.count(active_region) == 1) {
        throw HaplotypeOverflowError {active_region};
    } else {
        previous_holdout_regions_.emplace(active_region);
    }
    
    auto active_alleles = copy_contained(alleles_, active_region);
    
    assert(!active_alleles.empty());
    
    auto active_regions = extract_regions(active_alleles);
    
    active_regions.erase(std::unique(std::begin(active_regions), std::end(active_regions)),
                         std::end(active_regions));
    
    std::vector<std::pair<GenomicRegion, unsigned>> interaction_counts {};
    interaction_counts.reserve(active_regions.size());
    
    std::transform(std::cbegin(active_regions), std::cend(active_regions),
                   std::back_inserter(interaction_counts),
                   [&active_alleles] (const auto& region) {
                       return std::make_pair(region, count_overlapped(active_alleles, region));
                   });
    
    std::sort(std::begin(interaction_counts), std::end(interaction_counts),
              [] (const auto& lhs, const auto& rhs) { return lhs.second < rhs.second; });
    
    auto new_active_region = current_active_region_;
    
    do {
        const auto& holdout_region = interaction_counts.back().first;
        
        // it is safe to take base iterators as alleles with the same region must be adjacent
        auto contained = bases(contained_range(active_alleles, holdout_region));
        
        const auto is_holdout = [&holdout_region] (const auto& allele) { return is_same_region(allele, holdout_region); };
        
        const auto it1 = std::find_if(std::cbegin(contained), std::cend(contained), is_holdout);
        assert(it1 != std::cend(contained));
        const auto it2 = std::find_if_not(std::next(it1), std::cend(contained), is_holdout);
        
        holdout_set_.insert(it1, it2);
        active_alleles.erase(it1, it2);
        
        new_active_region = walker_.walk(head_region(current_active_region_), reads_, active_alleles);
        
        interaction_counts.pop_back();
    } while (!active_alleles.empty()
             && std::exp2(count_overlapped(active_alleles, new_active_region)) > haplotype_limits_.hard_max);
    
    if (DEBUG_MODE) {
        Logging::DebugLogger log {};
        stream(log) << "Adding " << holdout_set_.size() << " alleles to the holdout set";
    }
    
    current_holdout_region_ = encompassing_region(holdout_set_);
    
    alleles_.erase_all(std::cbegin(holdout_set_), std::cend(holdout_set_));
}

bool HaplotypeGenerator::try_reintroducing_holdout_set()
{
    if (holdout_set_.empty()) return false;
    
    if (ends_before(current_active_region_, *current_holdout_region_)
        && has_overlapped(alleles_, right_overhang_region(*current_holdout_region_, current_active_region_))) {
        return false;
    }
    
    if (DEBUG_MODE) {
        Logging::DebugLogger log {};
        stream(log) << "Reintroducing " << holdout_set_.size() << " holdout alleles";
    }
    
//    std::cout << "haplotypes before splice" << std::endl;
//    for (const auto& haplotype : tree_.extract_haplotypes()) {
//        ::debug::print_alleles(haplotype); std::cout << '\n';
//    }
    
    tree_.splice(Allele {"22", 49189922, "T"});
    
//    for (const auto& allele : holdout_set_) {
//        tree_.splice(allele);
//    }
    
//    std::cout << "haplotypes after splice" << std::endl;
//    for (const auto& haplotype : tree_.extract_haplotypes()) {
//        ::debug::print_alleles(haplotype); std::cout << '\n';
//    }
    
//    // We want to add the holdout alleles to the existing tree without regenerating already
//    // pruned haplotypes. HaplotypeTree does not support threading alleles into the existing structure
//    // (only extension), so we will need to re-construct a full tree (with the holdout alleles)
//    // and then prune haplotypes that are made up of only non-holdout alleles and
//    // were not in the current tree
//    
//    auto existing_haplotypes = tree_.extract_haplotypes();
//    
//    auto holdout_interaction_region = closed_region(*current_holdout_region_, current_active_region_);
//    
//    tree_.remove_overlapped(holdout_interaction_region);
//    
//    alleles_.insert(std::begin(holdout_set_), std::end(holdout_set_));
//    
//    auto holdout_interation_alleles = overlap_range(alleles_, holdout_interaction_region);
//    
//    extend_tree(holdout_interation_alleles, tree_);
//    
//    auto all_possible_haplotypes = tree_.extract_haplotypes();
//    
//    std::sort(std::begin(existing_haplotypes), std::end(existing_haplotypes));
//    std::sort(std::begin(all_possible_haplotypes), std::end(all_possible_haplotypes));
//    
//    // we only need to consider unqiue haplotypes as we will remove all duplicates from the tree
//    existing_haplotypes.erase(std::unique(std::begin(existing_haplotypes),
//                                          std::end(existing_haplotypes)),
//                              std::end(existing_haplotypes));
//    all_possible_haplotypes.erase(std::unique(std::begin(all_possible_haplotypes),
//                                              std::end(all_possible_haplotypes)),
//                                  std::end(all_possible_haplotypes));
//    
//    std::deque<Haplotype> haplotypes_to_remove {};
//    
//    std::set_difference(std::make_move_iterator(std::begin(all_possible_haplotypes)),
//                        std::make_move_iterator(std::end(all_possible_haplotypes)),
//                        std::cbegin(existing_haplotypes),
//                        std::cend(existing_haplotypes),
//                        std::back_inserter(haplotypes_to_remove));
//    
//    all_possible_haplotypes.clear();
//    all_possible_haplotypes.shrink_to_fit();
//    existing_haplotypes.clear();
//    existing_haplotypes.shrink_to_fit();
//    
//    auto it = std::remove_if(std::begin(haplotypes_to_remove), std::end(haplotypes_to_remove),
//                             [this] (const auto& haplotype) {
//                                 return std::any_of(std::cbegin(holdout_set_), std::cend(holdout_set_),
//                                                    [&haplotype] (const auto& allele) {
//                                                        return haplotype.contains_exact(allele);
//                                                    });
//                              });
//    
//    haplotypes_to_remove.erase(it, std::end(haplotypes_to_remove));
//    
//    prune_all(haplotypes_to_remove, tree_);
    
    holdout_set_.clear();
    holdout_set_.shrink_to_fit();
    
    current_holdout_region_ = boost::none;
    
    return true;
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
    const auto overlapped = overlap_range(alleles_, current_active_region_);
    
    // We want to keep haplotypes as small as possible, while allowing sufficient flanking
    // reference sequence for full read re-mapping and alignment (i.e. the read must be
    // contained by the haplotype). Note the sum of the indel sizes may not be sufficient
    // as the candidate generator may not propopse all variation in the original reads.
    const auto additional_padding = 2 * sum_indel_sizes(overlapped) + min_pad_;
    
    if (has_overlapped(reads_.get(), current_active_region_)) {
        const auto& lhs_read = *leftmost_overlapped(reads_.get(), current_active_region_);
        const auto& rhs_read = *rightmost_overlapped(reads_.get(), current_active_region_);
        
        const auto unpadded_region = encompassing_region(lhs_read, rhs_read);
        
        if (mapped_begin(lhs_read) < additional_padding / 2) {
            const auto lhs_padding = mapped_begin(lhs_read);
            const auto rhs_padding = additional_padding - lhs_padding;
            return expand(unpadded_region, lhs_padding, rhs_padding);
        }
        
        return expand(unpadded_region, additional_padding / 2);
    }
    
    return expand(current_active_region_, additional_padding / 2);
}

// Builder

HaplotypeGenerator::Builder& HaplotypeGenerator::Builder::set_soft_max_haplotypes(unsigned n) noexcept
{
    haplotype_limits_.soft_max = n;
    haplotype_limits_.holdout_max = std::max(haplotype_limits_.holdout_max, n + 1);
    haplotype_limits_.hard_max = std::max(haplotype_limits_.hard_max, n + 1);
    return *this;
}

HaplotypeGenerator::Builder& HaplotypeGenerator::Builder::set_hard_max_haplotypes(unsigned n) noexcept
{
    haplotype_limits_.hard_max = n;
    return *this;
}

HaplotypeGenerator::Builder& HaplotypeGenerator::Builder::set_lagging_policy(LaggingPolicy policy) noexcept
{
    lagging_policy_ = policy;
    return *this;
}

HaplotypeGenerator::Builder& HaplotypeGenerator::Builder::set_min_pad(Haplotype::SizeType val) noexcept
{
    min_pad_ = val;
    return *this;
}

HaplotypeGenerator HaplotypeGenerator::Builder::build(const ReferenceGenome& reference, const GenomicRegion& window,
                         const MappableFlatSet<Variant>& candidates,
                         const ReadMap& reads) const
{
    return HaplotypeGenerator {
        window, reference, candidates, reads, haplotype_limits_, lagging_policy_, min_pad_
    };
}
} // namespace Octopus
