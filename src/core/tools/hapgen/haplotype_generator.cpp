// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "haplotype_generator.hpp"

#include <algorithm>
#include <deque>
#include <iterator>
#include <numeric>
#include <cmath>
#include <cassert>

#include <boost/range/iterator_range.hpp>

#include "core/types/variant.hpp"
#include "core/types/haplotype.hpp"
#include "concepts/mappable.hpp"
#include "utils/mappable_algorithms.hpp"
#include "utils/append.hpp"

#include <iostream> // DEBUG
#include "timers.hpp"

#define _unused(x) ((void)(x))

namespace octopus { namespace coretools {

// HaplotypeOverflow

HaplotypeGenerator::HaplotypeOverflow::HaplotypeOverflow(GenomicRegion region, std::size_t size)
: runtime_error {"HaplotypeOverflowError"}
, region_ {std::move(region)}
, size_ {size}
, message_ {}
{}

const char* HaplotypeGenerator::HaplotypeOverflow::what() const noexcept
{
    return runtime_error::what();
}

const GenomicRegion& HaplotypeGenerator::HaplotypeOverflow::region() const noexcept
{
    return region_;
}

std::size_t HaplotypeGenerator::HaplotypeOverflow::size() const noexcept
{
    return size_;
}

// HaplotypeGenerator

namespace {

auto max_included(const unsigned max_haplotypes) noexcept
{
    return 2 * static_cast<unsigned>(std::max(1.0, std::log2(std::max(max_haplotypes, 1u)))) - 1;
}

auto get_walker_policy(const HaplotypeGenerator::Policies::Extension policy) noexcept
{
    using HGP = HaplotypeGenerator::Policies::Extension;
    using GWP = GenomeWalker::ExtensionPolicy;
    switch (policy) {
        case HGP::conservative: return GWP::includeIfWithinReadLengthOfFirstIncluded;
        case HGP::normal: return GWP::includeIfAllSamplesSharedWithFrontier;
        case HGP::optimistic: return GWP::includeIfAnySampleSharedWithFrontier;
        case HGP::aggressive: return GWP::noLimit;
        default: return GWP::includeIfAllSamplesSharedWithFrontier; // prevents compiler warning
    }
}

auto get_walker_policy(const HaplotypeGenerator::Policies::Lagging policy) noexcept
{
    using HGP = HaplotypeGenerator::Policies::Lagging;
    using GWP = GenomeWalker::IndicatorPolicy;
    switch (policy) {
        case HGP::none: return GWP::includeNone;
        case HGP::conservative:
        case HGP::moderate:
        case HGP::normal: return GWP::includeIfSharedWithNovelRegion;
        case HGP::aggressive: return GWP::includeIfLinkableToNovelRegion;
        default: return GWP::includeIfSharedWithNovelRegion; // prevents compiler warning
    }
}

auto decompose(const MappableFlatSet<Variant>& variants)
{
    std::vector<Allele> alleles {};
    alleles.reserve(2 * variants.size());
    for (const auto& variant : variants) {
        alleles.push_back(variant.ref_allele());
        alleles.push_back(variant.alt_allele());
    }
    return MappableFlatSet<Allele> {
        std::make_move_iterator(std::begin(alleles)),
        std::make_move_iterator(std::end(alleles))
    };
}

}

namespace debug {
    template <typename Range>
    void print_active_alleles(const Range& alleles, const GenomicRegion& active_region);
} // namespace debug

// public members

HaplotypeGenerator::HaplotypeGenerator(const ReferenceGenome& reference, const MappableFlatSet<Variant>& candidates,
                                       const ReadMap& reads, Policies policies)
try
: policies_ {std::move(policies)}
, tree_ {contig_name(candidates.front()), reference}
, default_walker_ {
    max_included(policies_.haplotype_limits.target),
    GenomeWalker::IndicatorPolicy::includeNone,
    get_walker_policy(policies.extension)
}
, holdout_walker_ {
    max_included(policies_.haplotype_limits.target),
    GenomeWalker::IndicatorPolicy::includeAll,
    get_walker_policy(policies.extension)
  }
, lagged_walker_ {}
, alleles_ {decompose(candidates)}
, reads_ {reads}
, next_active_region_ {}
, active_holdouts_ {}
, holdout_region_ {}
, previous_holdout_regions_ {}
, debug_log_ {logging::get_debug_log()}
, trace_log_ {logging::get_trace_log()}
{
    rightmost_allele_ = alleles_.rightmost();
    active_region_ = head_region(alleles_.leftmost());
    if (active_region_.begin() != 0) {
        active_region_ = shift(active_region_, -1);
    }
    if (policies.lagging != Policies::Lagging::none) {
        lagged_walker_ = GenomeWalker {
            max_included(policies_.haplotype_limits.target),
            get_walker_policy(policies_.lagging),
            get_walker_policy(policies_.extension)
        };
    }
}
catch (...) {
    if (candidates.empty()) {
        throw std::runtime_error {"HaplotypeGenerator: not supplied with any candidates"};
    }
    throw;
}

namespace {

template <typename ForwardIt, typename MappableTp>
auto find_first_exact_overlap(ForwardIt first, ForwardIt last, const MappableTp& mappable)
{
    return std::find_if(first, last, [&] (const auto& m) { return is_same_region(m, mappable); });
}

template <typename Container, typename MappableTp>
auto find_first_exact_overlap(const Container& container, const MappableTp& mappable)
{
    return find_first_exact_overlap(std::cbegin(container), std::cend(container), mappable);
}

bool has_rhs_sandwich_insertion(const MappableFlatSet<Allele>& alleles, const GenomicRegion& active_region) {
    const auto insertion_region = tail_region(active_region);
    const auto active_tails = overlap_range(alleles, insertion_region);
    if (empty(active_tails)) return false;
    const auto insertion_itr = find_first_exact_overlap(active_tails, insertion_region);
    if (insertion_itr != std::cend(active_tails)) {
        return begins_before(active_tails.front(), insertion_region) && ends_before(insertion_region, active_tails.back());
    } else {
        return false;
    }
}

} // namespace

HaplotypeGenerator::HaplotypePacket HaplotypeGenerator::generate()
{
    if (alleles_.empty()) return std::make_tuple(std::vector<Haplotype> {}, active_region_, boost::none);
    populate_tree();
    const auto haplotype_region = calculate_haplotype_region();
    assert(contains(haplotype_region, active_region_));
    auto haplotypes = tree_.extract_haplotypes(haplotype_region);
    if (!(is_lagging_enabled() || in_holdout_mode())) tree_.clear();
    return std::make_tuple(std::move(haplotypes), active_region_, holdout_region_);
}

boost::optional<GenomicRegion> HaplotypeGenerator::peek_next_active_region() const
{
    if (in_holdout_mode()) return boost::none;
    update_next_active_region();
    return *next_active_region_;
}

void HaplotypeGenerator::clear_progress() noexcept
{
    tree_.clear();
    reset_next_active_region();
    if (in_holdout_mode()) {
        clear_holdouts();
    }
}

void HaplotypeGenerator::jump(GenomicRegion region)
{
    clear_progress();
    progress(std::move(region));
}

bool HaplotypeGenerator::removal_has_impact() const
{
    if (in_holdout_mode()) return true;
    if (!is_lagging_enabled() || contains(active_region_, rightmost_allele_)) return false;
    const auto max_lagged_region = lagged_walker_->walk(active_region_, reads_, alleles_);
    return overlaps(max_lagged_region, active_region_);
}

unsigned HaplotypeGenerator::max_removal_impact() const
{
    if (in_holdout_mode()) return tree_.num_haplotypes();
    if (!is_lagging_enabled() || contains(active_region_, rightmost_allele_)) return 0;
    const auto max_lagged_region = lagged_walker_->walk(active_region_, reads_, alleles_);
    if (!overlaps(max_lagged_region, active_region_)) return 0;
    const auto novel_region = right_overhang_region(max_lagged_region, active_region_);
    const auto num_novel_alleles = count_overlapped(alleles_, novel_region);
    if (num_novel_alleles == 0) return 0;
    const auto max_new_haplotypes = std::max(static_cast<unsigned>(std::exp2(num_novel_alleles / 2)), 1u);
    const auto num_leftover_haplotypes = policies_.haplotype_limits.target / max_new_haplotypes;
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
    if (in_holdout_mode()) return true;
    if (!is_lagging_enabled()) return false;
    const auto next_lagged_region = lagged_walker_->walk(active_region_, reads_, alleles_);
    return overlaps(active_region_, next_lagged_region);
}

void HaplotypeGenerator::reset_next_active_region() const noexcept
{
    next_active_region_ = boost::none;
}

void HaplotypeGenerator::update_next_active_region() const
{
    if (!next_active_region_) {
        if (is_lagging_enabled() || in_holdout_mode()) {
            // If we are in holdout mode then lagging is required
            update_lagged_next_active_region();
        } else {
            next_active_region_ = default_walker_.walk(active_region_, reads_, alleles_);
        }
    }
    assert(in_holdout_mode() || active_region_ <= *next_active_region_);
}

namespace {

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
    const auto last = crend(passed_alleles);
    const auto& last_passed_allele = passed_alleles.back();
    const auto it = std::find_if_not(std::next(crbegin(passed_alleles)), last,
                                     [&last_passed_allele] (const auto& allele) {
                                         return is_same_region(allele, last_passed_allele);
                                     });
    return it != last && is_position(*it);
}

template <typename T>
void pop_front(std::vector<T>& v) {
    assert(!v.empty());
    v.erase(std::cbegin(v));
}

} // namespace

GenomicRegion HaplotypeGenerator::find_max_lagged_region() const
{
    if (in_holdout_mode()) {
        return holdout_walker_.walk(active_region_, reads_, alleles_);
    } else {
        return lagged_walker_->walk(active_region_, reads_, alleles_);
    }
}

namespace {

bool try_extend_tree_without_removal(HaplotypeTree& tree, const GenomicRegion& active_region,
                                     const MappableFlatSet<Allele>& alleles,
                                     const GenomicRegion& max_lagged_region,
                                     const unsigned haplotype_limit)
{
    const auto novel_region = right_overhang_region(max_lagged_region, active_region);
    const auto novel_alleles = overlap_range(alleles, novel_region);
    const auto last_added = extend_tree_until(novel_alleles, tree, haplotype_limit);
    if (last_added == std::cend(novel_alleles)) {
        return true;
    } else {
        tree.clear(novel_region); // undo previous extension
        return false;
    }
}

void safe_clear_passed(HaplotypeTree& tree, const GenomicRegion& active_region,
                       const MappableFlatSet<Allele>& alleles, const GenomicRegion& max_lagged_region)
{
    const auto passed_region = left_overhang_region(active_region, max_lagged_region);
    const auto passed_alleles = overlap_range(alleles, passed_region);
    if (can_remove_entire_passed_region(active_region, max_lagged_region, passed_alleles)) {
        tree.clear(passed_region);
    } else if (requires_staged_removal(passed_alleles)) {
        const auto first_removal_region = expand_rhs(passed_region, -1);
        tree.clear(first_removal_region);
        tree.clear(tail_region(first_removal_region));
    } else {
        tree.clear(expand_rhs(passed_region, -1));
    }
}

void prune_indicators(HaplotypeTree& tree, std::vector<GenomicRegion>& indicator_regions,
                      const std::size_t target_tree_size)
{
    auto itr = std::find_if(std::cbegin(indicator_regions), std::cend(indicator_regions),
                            [&tree, target_tree_size](const auto& region) {
                                if (tree.num_haplotypes() > target_tree_size) {
                                    tree.clear(region);
                                    return false;
                                } else {
                                    return true;
                                }
                            });
    indicator_regions.erase(std::cbegin(indicator_regions), itr);
}

template <typename Range>
unsigned extend_novel(HaplotypeTree& tree, const std::vector<GenomicRegion>& novel_regions,
                      const Range& novel_alleles, const HaplotypeGenerator::Policies::HaplotypeLimits& limits)
{
    unsigned num_novel_regions_added {0};
    for (const auto& region : novel_regions) {
        const auto interacting_alleles = contained_range(novel_alleles, region);
        const auto last_added = extend_tree_until(interacting_alleles, tree, limits.overflow);
        if (last_added != std::cend(interacting_alleles)) {
            break;
        }
        ++num_novel_regions_added;
        if (tree.num_haplotypes() > limits.target) {
            if (num_novel_regions_added > 1
                && !(num_novel_regions_added == novel_regions.size()
                     && tree.num_haplotypes() < limits.holdout)) {
                tree.clear(region);
                --num_novel_regions_added;
            }
            break;
        } else if (tree.num_haplotypes() == limits.target) {
            break;
        }
    }
    return num_novel_regions_added;
}

bool has_duplicate_insertion(const std::vector<GenomicRegion>& indicator_regions,
                             const std::vector<GenomicRegion>& novel_regions)
{
    return is_empty(novel_regions.front()) && are_adjacent(indicator_regions.back(), novel_regions.front());
}

void remove_duplicate_novels(const std::vector<GenomicRegion>& indicator_regions,
                             std::vector<GenomicRegion>& novel_regions)
{
    if (has_duplicate_insertion(indicator_regions, novel_regions)) {
        pop_front(novel_regions);
    }
}

std::size_t get_reduction_factor(HaplotypeGenerator::Policies::Lagging policy) noexcept
{
    using HGP = HaplotypeGenerator::Policies::Lagging;
    switch (policy) {
        case HGP::none:
        case HGP::conservative: return 0;
        case HGP::moderate: return 2;
        case HGP::normal: return 3;
        case HGP::aggressive: return 4;
        default: return 3; // prevents compiler warning
    }
}

std::size_t get_target_tree_size(const std::size_t curr_tree_size, const std::size_t target_tree_size,
                                 const std::size_t num_mutually_exclusive_indicator_regions,
                                 const std::size_t num_mutually_exclusive_novel_regions,
                                 const std::size_t reduction_factor)
{
    if (curr_tree_size == 0) return 0;
    const auto effective_log_space = static_cast<std::size_t>(std::log2(std::max(target_tree_size / curr_tree_size, std::size_t {1})));
    if (effective_log_space >= num_mutually_exclusive_novel_regions) {
        return curr_tree_size;
    } else {
        const auto exponent = std::min(num_mutually_exclusive_indicator_regions,
                                        (std::max(num_mutually_exclusive_novel_regions, reduction_factor)
                                         / std::max(reduction_factor, std::size_t {1})) - 1);
        return std::min(static_cast<std::size_t>(target_tree_size / std::pow(2, exponent)), curr_tree_size);
    }
}

std::size_t get_target_tree_size(const std::size_t curr_tree_size, const std::size_t target_tree_size,
                                 const std::size_t num_mutually_exclusive_indicator_regions,
                                 const std::size_t num_mutually_exclusive_novel_regions,
                                 HaplotypeGenerator::Policies::Lagging policy)
{
    const auto reduction_factor = get_reduction_factor(policy);
    return get_target_tree_size(curr_tree_size, target_tree_size, num_mutually_exclusive_indicator_regions,
                                num_mutually_exclusive_novel_regions, reduction_factor);
}

auto max_ref_distance(const GenomicRegion& passed_region, const HaplotypeTree& tree)
{
    HaplotypeTree::HaplotypeLength min_length, max_length;
    std::tie(min_length, max_length) = minmax_haplotype_lengths(tree, passed_region);
    if (min_length == max_length) {
        return HaplotypeTree::HaplotypeLength {0};
    } else {
        if (min_length <= size(passed_region) && size(passed_region) <= max_length) {
            return std::max(max_length - size(passed_region), size(passed_region) - min_length);
        } else if (min_length <= size(passed_region)) {
            return size(passed_region) - min_length;
        } else if (size(passed_region) <= max_length) {
            return max_length - size(passed_region);
        } else {
            return HaplotypeTree::HaplotypeLength {0};
        }
    }
}

bool requires_exclusion_zone(const GenomicRegion& passed_region, const HaplotypeTree& tree,
                             const Haplotype::MappingDomain::Size min_flank_pad)
{
    return max_ref_distance(passed_region, tree) >= min_flank_pad;
}

boost::optional<GenomicRegion>
get_exclusion_zone(const GenomicRegion& passed_region, const HaplotypeTree& tree,
                   const Haplotype::MappingDomain::Size min_flank_pad,
                   const ReadMap& reads)
{
    if (requires_exclusion_zone(passed_region, tree, min_flank_pad)) {
        return encompassing_overlap_region(reads, passed_region);
    } else {
        return boost::none;
    }
}

} // namespace

void HaplotypeGenerator::update_lagged_next_active_region() const
{
    if (contains(active_region_, rightmost_allele_)) {
        next_active_region_ = shift(tail_region(rightmost_allele_), 2); // Nothing more to do
        return;
    }
    const auto max_lagged_region = find_max_lagged_region();
    assert(has_contained(alleles_, max_lagged_region));
    if (!overlaps(active_region_, max_lagged_region)) {
        next_active_region_ = std::move(max_lagged_region);
    } else {
        HaplotypeTree test_tree {tree_}; // use a temporary tree to see how much we can lag
        if (begins_before(active_region_, max_lagged_region)) {
            if (try_extend_tree_without_removal(test_tree, active_region_, alleles_,  max_lagged_region,
                                                policies_.haplotype_limits.target)) {
                next_active_region_ = test_tree.encompassing_region();
                return;
            } else if (!in_holdout_mode()) {
                safe_clear_passed(test_tree, active_region_, alleles_,  max_lagged_region);
            }
        }
        // overlap_range is required for novel alleles as any holdout alleles
        // just reintroduced may overlap with the indicator and novel region.
        const auto novel_region  = right_overhang_region(max_lagged_region, active_region_);
        const auto novel_alleles = overlap_range(alleles_, novel_region);
        assert(!novel_alleles.empty());
        auto mutually_exclusive_novel_regions = extract_mutually_exclusive_regions(novel_alleles);
        // contained_range is required for indicators as previous active alleles must be contained
        // within the reported active region.
        const auto indicator_region  = *overlapped_region(active_region_, max_lagged_region);
        const auto indicator_alleles = contained_range(alleles_, indicator_region);
        if (!indicator_alleles.empty()) {
            auto mutually_exclusive_indicator_regions = extract_mutually_exclusive_regions(indicator_alleles);
            // Although HaplotypeTree can handle duplicate extension, if the tree is not in a complete
            // state (i.e. it has been pruned) then adding duplicates will corrupt the tree.
            remove_duplicate_novels(mutually_exclusive_indicator_regions, mutually_exclusive_novel_regions);
            if (!in_holdout_mode()) {
                const auto novel_overlap_region = encompassing_region(novel_alleles);
                while (!mutually_exclusive_indicator_regions.empty()
                       && overlaps(mutually_exclusive_indicator_regions.back(), novel_overlap_region)) {
                    mutually_exclusive_indicator_regions.pop_back();
                }
                auto target_tree_size = get_target_tree_size(test_tree.num_haplotypes(),
                                                             policies_.haplotype_limits.target,
                                                             mutually_exclusive_indicator_regions.size(),
                                                             mutually_exclusive_novel_regions.size(),
                                                             policies_.lagging);
                prune_indicators(test_tree, mutually_exclusive_indicator_regions, target_tree_size);
                if (!mutually_exclusive_indicator_regions.empty()) {
                    const auto new_passed_region = left_overhang_region(active_region_, mutually_exclusive_indicator_regions.front());
                    auto passed_alleles = contained_range(alleles_, new_passed_region);
                    while (!empty(passed_alleles) && overlaps(passed_alleles.back(), mutually_exclusive_indicator_regions.front())) {
                        passed_alleles.advance_end(-1);
                    }
                    if (!empty(passed_alleles)) {
                        const auto passed_allele_region = encompassing_region(passed_alleles);
                        const auto exclusion_zone = get_exclusion_zone(passed_allele_region, tree_,
                                                                       policies_.exclusion_threshold, reads_);
                        if (exclusion_zone) {
                            if (debug_log_) {
                                const auto new_indicator_region = encompassing_region(mutually_exclusive_indicator_regions);
                                stream(*debug_log_) << "Encountered indicator exclusion zone " << *exclusion_zone
                                                    << " while trying to use indicator region " << new_indicator_region;
                            }
                            test_tree.clear(*exclusion_zone);
                        }
                    }
                }
            }
        }
        if (in_holdout_mode() && ends_before(top_holdout_region(), novel_region)) {
            const auto holdout_novel_region = closed_region(novel_region, top_holdout_region());
            mutually_exclusive_novel_regions = copy_overlapped(mutually_exclusive_novel_regions, holdout_novel_region);
        }
        const auto num_novel_regions_added = extend_novel(test_tree, mutually_exclusive_novel_regions,
                                                          novel_alleles, policies_.haplotype_limits);
        if (num_novel_regions_added > 0) {
            next_active_region_ = test_tree.encompassing_region();
        } else {
            assert(!mutually_exclusive_novel_regions.empty());
            if (test_tree.is_empty()) {
                next_active_region_ = mutually_exclusive_novel_regions.front();
            } else {
                next_active_region_ = closed_region(test_tree.encompassing_region(), mutually_exclusive_novel_regions.front());
            }
        }
        if (*next_active_region_ == active_region_) {
            next_active_region_ = default_walker_.walk(active_region_, reads_, alleles_);
        }
    }
}

void HaplotypeGenerator::progress(GenomicRegion to)
{
    if (to == active_region_) return;
    next_active_region_ = std::move(to);
    if (!in_holdout_mode()) {
        if (begins_before(active_region_, *next_active_region_)) {
            auto passed_region = left_overhang_region(active_region_, *next_active_region_);
            const auto passed_alleles = overlap_range(alleles_, passed_region);
            if (passed_alleles.empty()) return;
            if (can_remove_entire_passed_region(active_region_, *next_active_region_, passed_alleles)) {
                alleles_.erase_overlapped(passed_region);
                tree_.clear(passed_region);
            } else if (requires_staged_removal(passed_alleles)) {
                // We need to be careful here as insertions adjacent to passed_region are
                // considered overlapped and would be wrongly erased if we erased the whole
                // region. But, we also want to clear all single base alleles left adjacent with
                // next_active_region_, as they have truly been passed.
                
                // This will erase everything to the left of the adjacent insertion, other than
                // the single base alleles adjacent with next_active_region_.
                const auto first_removal_region = expand_rhs(passed_region, -1);
                alleles_.erase_overlapped(first_removal_region);
                // This will erase the remaining single base alleles in passed_region, but not the
                // insertions in next_active_region_.
                const auto second_removal_region = tail_region(first_removal_region);
                alleles_.erase_overlapped(second_removal_region);
                
                if (is_after(*next_active_region_, active_region_)) {
                    assert(tree_.is_empty() || contains(active_region_, tree_.encompassing_region()));
                    tree_.clear();
                } else {
                    tree_.clear(first_removal_region);
                    tree_.clear(second_removal_region);
                }
            } else {
                const auto removal_region = expand_rhs(passed_region, -1);
                alleles_.erase_overlapped(removal_region);
                tree_.clear(removal_region);
            }
            if (overlaps(passed_region, rightmost_allele_) && !alleles_.empty()) {
                rightmost_allele_ = alleles_.rightmost();
            }
        }
    }
}

void HaplotypeGenerator::populate_tree()
{
    if (in_holdout_mode() && can_reintroduce_holdouts()) {
        populate_tree_with_holdouts();
    } else {
        populate_tree_with_novel_alleles();
    }
}

void HaplotypeGenerator::populate_tree_with_novel_alleles()
{
    update_next_active_region();
    if (is_after(*next_active_region_, rightmost_allele_)) {
        // We are done
        tree_.clear();
        active_region_ = *next_active_region_;
        return;
    }
    progress(*next_active_region_);
    auto novel_active_region = *next_active_region_;
    if (!tree_.is_empty()) {
        novel_active_region = right_overhang_region(*next_active_region_, active_region_);
    }
    auto novel_active_alleles = overlap_range(alleles_, novel_active_region);
    auto last_added_novel_itr = extend_tree_until(novel_active_alleles, tree_, policies_.haplotype_limits.holdout);
    if (last_added_novel_itr != std::cend(novel_active_alleles)) {
        reset_next_active_region();
        const auto active_region_before_holdout = active_region_;
        auto next_holdout_region = novel_active_region;
        while (last_added_novel_itr != std::cend(novel_active_alleles) && can_try_extracting_holdouts(next_holdout_region)) {
            if (!try_extract_holdouts(next_holdout_region)) {
                break;
            }
            tree_.clear(novel_active_region);
            active_region_ = active_region_before_holdout;
            update_next_active_region();
            active_region_ = *std::move(next_active_region_);
            reset_next_active_region();
            novel_active_region = right_overhang_region(active_region_, active_region_before_holdout);
            novel_active_alleles = overlap_range(alleles_, novel_active_region);
            last_added_novel_itr = extend_tree_until(novel_active_alleles, tree_, policies_.haplotype_limits.overflow);
            if (overlaps(active_region_, top_holdout_region())) {
                next_holdout_region = *overlapped_region(active_region_, top_holdout_region());
            } else {
                next_holdout_region = novel_active_region;
            }
        }
        if (last_added_novel_itr != std::cend(novel_active_alleles)) {
            last_added_novel_itr = extend_tree_until(last_added_novel_itr, std::cend(novel_active_alleles), tree_,
                                                     policies_.haplotype_limits.overflow);
            active_region_ = tree_.encompassing_region();
            if (last_added_novel_itr != std::cend(novel_active_alleles)) {
                throw HaplotypeOverflow {active_region_, tree_.num_haplotypes()};
            }
        }
    } else {
        active_region_ = *std::move(next_active_region_);
        reset_next_active_region();
    }
}

void HaplotypeGenerator::populate_tree_with_holdouts()
{
    reintroduce_holdouts();
    if (tree_.num_haplotypes() > policies_.haplotype_limits.overflow) {
        throw HaplotypeOverflow {active_region_, tree_.num_haplotypes()};
    }
    active_region_ = tree_.encompassing_region();
    if (!in_holdout_mode() && has_rhs_sandwich_insertion(alleles_, active_region_)) {
        resolve_sandwich_insertion();
    }
    reset_next_active_region();
}

bool HaplotypeGenerator::in_holdout_mode() const noexcept
{
    assert((holdout_region_ && !active_holdouts_.empty()) || ((!holdout_region_ && active_holdouts_.empty())));
    return !active_holdouts_.empty();
}

const GenomicRegion& HaplotypeGenerator::top_holdout_region() const
{
    assert(in_holdout_mode());
    return active_holdouts_.top().region;
}

bool HaplotypeGenerator::can_try_extracting_holdouts(const GenomicRegion& region) const noexcept
{
    return active_holdouts_.size() < policies_.max_holdout_depth && has_overlapped(alleles_, region);
}

namespace debug {

template <typename S, typename Contianer>
void print_new_holdouts(S&& stream, const Contianer& alleles)
{
    stream << "Adding " << alleles.size() << " alleles to the holdout stack:" << '\n';
    for (const auto& allele : alleles) {
        stream << allele << '\n';
    }
}

} // namespace debug

template <typename Range>
unsigned estimate_num_haplotype(const Range& alleles, const GenomicRegion& next_active_region)
{
    return std::exp2(count_overlapped(alleles, next_active_region));
}

template <typename Range>
bool require_more_holdouts(const Range& alleles, const GenomicRegion& next_active_region,
                           const unsigned haplotype_limit)
{
    return !alleles.empty() && estimate_num_haplotype(alleles, next_active_region) > haplotype_limit;
}

namespace {

auto extract_unique_regions(const MappableFlatSet<Allele>& alleles)
{
    auto result = extract_regions(alleles);
    result.erase(std::unique(std::begin(result), std::end(result)), std::end(result));
    return result;
}

struct InteractionCount
{
    GenomicRegion region;
    std::size_t count;
};

auto get_interaction_counts(const MappableFlatSet<Allele>& alleles)
{
    const auto allele_regions = extract_unique_regions(alleles);
    std::vector<InteractionCount> result {};
    result.reserve(result.size());
    std::transform(std::cbegin(allele_regions), std::cend(allele_regions), std::back_inserter(result),
                   [&alleles] (const auto& region) -> InteractionCount {
                       return {region, count_overlapped(alleles, region)};
                   });
    std::sort(std::begin(result), std::end(result),
              [] (const auto& lhs, const auto& rhs) {
                  if (lhs.count == rhs.count) {
                      // Prefer holding out larger regions
                      return size(lhs.region) > size(rhs.region);
                  }
                  return lhs.count > rhs.count;
              });
    return result;
}

using HoldoutRange = boost::iterator_range<MappableFlatSet<Allele>::const_iterator>;

HoldoutRange get_holdout_range(const MappableFlatSet<Allele>& alleles, const GenomicRegion& holdout_region)
{
    // It is safe to take base iterators as alleles with the same region must be adjacent
    const auto contained = bases(contained_range(alleles, holdout_region));
    const auto is_holdout = [&holdout_region] (const auto& allele) {
        return is_same_region(allele, holdout_region);
    };
    const auto first_holdout = std::find_if(std::cbegin(contained), std::cend(contained), is_holdout);
    assert(first_holdout != std::cend(contained));
    const auto last_holdout = std::find_if_not(std::next(first_holdout), std::cend(contained), is_holdout);
    return boost::make_iterator_range(first_holdout, last_holdout);
}

void remove_move(MappableFlatSet<Allele>& src, const HoldoutRange& holdouts, std::deque<Allele>& dst)
{
    dst.insert(std::cend(dst), std::cbegin(holdouts), std::cend(holdouts));
    src.erase(std::cbegin(holdouts), std::cend(holdouts));
}

auto demote_all(const std::vector<GenomicRegion>& regions)
{
    std::vector<ContigRegion> result(regions.size());
    std::transform(std::cbegin(regions), std::cend(regions), std::begin(result),
                   [] (const auto& region) { return region.contig_region(); });
    return result;
}

} // namespace

namespace debug {

void log_new_holdouts(const std::deque<Allele>& holdouts, boost::optional<logging::DebugLogger>& log)
{
    if (log) {
        debug::print_new_holdouts(stream(*log), holdouts);
    }
}

} // namespace debug

bool HaplotypeGenerator::try_extract_holdouts(GenomicRegion region)
{
    auto new_holdout_set = propose_new_holdout_set(std::move(region));
    if (new_holdout_set) {
        std::deque<Allele> new_holdout_alleles {};
        auto new_holdout_regions = demote_all(extract_mutually_exclusive_regions(new_holdout_set->alleles));
        if (previous_holdout_regions_.count(new_holdout_regions) == 0) {
            utils::append(new_holdout_set->alleles, new_holdout_alleles);
            active_holdouts_.push(std::move(*new_holdout_set));
            previous_holdout_regions_.insert(std::move(new_holdout_regions));
        }
        if (!new_holdout_alleles.empty()) {
            std::sort(std::begin(new_holdout_alleles), std::end(new_holdout_alleles));
            debug::log_new_holdouts(new_holdout_alleles, debug_log_);
            auto new_holdout_region = encompassing_region(new_holdout_alleles);
            if (holdout_region_) {
                holdout_region_ = encompassing_region(*holdout_region_, new_holdout_region);
            } else {
                holdout_region_ = std::move(new_holdout_region);
            }
            alleles_.erase_all(std::cbegin(new_holdout_alleles), std::cend(new_holdout_alleles));
            return true;
        } else {
            return false;
        }
    } else {
        return false;
    }
}

boost::optional<HaplotypeGenerator::HoldoutSet> HaplotypeGenerator::propose_new_holdout_set(GenomicRegion region) const
{
    assert(can_try_extracting_holdouts(region));
    auto active_alleles = copy_overlapped(alleles_, region);
    assert(!active_alleles.empty());
    auto interaction_counts = get_interaction_counts(active_alleles);
    const auto holdout_itr = std::find_if(std::cbegin(interaction_counts), std::cend(interaction_counts),
                                          [] (const auto& count) { return size(count.region) > 1; });
    if (holdout_itr != std::cend(interaction_counts)) {
        const auto holdouts = get_holdout_range(active_alleles, holdout_itr->region);
        assert(!holdouts.empty());
        return HoldoutSet {std::cbegin(holdouts), std::cend(holdouts), holdout_itr->region};
    } else {
        return boost::none;
    }
}

namespace {

template <typename Container, typename Mappable>
auto count_fully_contained(const Container& mappables, const Mappable& mappable)
{
    auto contained = contained_range(mappables, mappable);
    while (!empty(contained) && mapped_end(contained.front()) == mapped_begin(mappable)) {
        contained.advance_begin(1);
    }
    while (!empty(contained) && mapped_begin(contained.back()) == mapped_end(mappable)) {
        contained.advance_end(-1);
    }
    return size(contained);
}

template <typename Container, typename Mappable>
bool has_fully_contained(const Container& mappables, const Mappable& mappable)
{
    return count_fully_contained(mappables, mappable) > 0;
}

} // namespace

bool HaplotypeGenerator::can_reintroduce_holdouts() const noexcept
{
    assert(in_holdout_mode());
    if (!ends_before(active_region_, top_holdout_region())) {
        return true;
    } else {
        const auto remaining_holdout_region = right_overhang_region(top_holdout_region(), active_region_);
        return !has_fully_contained(alleles_, remaining_holdout_region);
    }
}

namespace debug {

template <typename S, typename Contianer>
void print_old_holdouts(S&& stream, const Contianer& alleles)
{
    stream << "Reintroducing " << alleles.size() << " holdout alleles:" << '\n';
    for (const auto& allele : alleles) {
        stream << allele << '\n';
    }
}

} // namespace debug

void HaplotypeGenerator::reintroduce_holdouts()
{
    assert(!active_holdouts_.empty());
    if (debug_log_) debug::print_old_holdouts(stream(*debug_log_), active_holdouts_.top().alleles);
    splice(active_holdouts_.top().alleles, tree_);
    if (ends_before(*holdout_region_, active_region_)) {
        auto extended_region = right_overhang_region(active_region_, *holdout_region_);
        extend_tree(contained_range(alleles_, extended_region), tree_);
    }
    alleles_.insert(std::make_move_iterator(std::begin(active_holdouts_.top().alleles)),
                    std::make_move_iterator(std::end(active_holdouts_.top().alleles)));
    active_holdouts_.pop();
    if (active_holdouts_.empty()) {
        holdout_region_ = boost::none;
    }
}

void HaplotypeGenerator::clear_holdouts() noexcept
{
    active_holdouts_ = decltype(active_holdouts_) {};
    holdout_region_ = boost::none;
}

void HaplotypeGenerator::resolve_sandwich_insertion()
{
    const auto active_alleles = bases(contained_range(alleles_, active_region_));
    auto first_required = std::cend(active_alleles);
    const auto last_required = find_next_mutually_exclusive(first_required, std::cend(alleles_));
    first_required = extend_tree_until(first_required, last_required, tree_, policies_.haplotype_limits.overflow);
    if (first_required != last_required) {
        throw HaplotypeOverflow {active_region_, tree_.num_haplotypes()};
    }
    active_region_ = tree_.encompassing_region();
}

template <typename Range>
auto sum_indel_sizes(const Range& alleles)
{
    return std::accumulate(std::cbegin(alleles), std::cend(alleles), std::size_t {0},
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
    const auto overlapped = overlap_range(alleles_, active_region_);
    // We want to keep haplotypes as small as possible, while allowing sufficient flanking
    // reference sequence for full read re-mapping and alignment (i.e. the read must be
    // contained by the haplotype). Note the sum of the indel sizes may not be sufficient
    // as the candidate generator may not propopse all variation in the original reads.
    const auto min_flank_padding = sum_indel_sizes(overlapped) + policies_.min_flank_pad;
    if (has_overlapped(reads_.get(), active_region_)) {
        const auto& lhs_read = *leftmost_overlapped(reads_.get(), active_region_);
        const auto& rhs_read = *rightmost_overlapped(reads_.get(), active_region_);
        const auto read_region = closed_region(lhs_read, rhs_read);
        GenomicRegion::Size lhs_expansion {}, rhs_expansion {};
        
        if (begins_before(read_region, active_region_)) {
            lhs_expansion = begin_distance(read_region, active_region_) + min_flank_padding;
        } else {
            const auto diff = static_cast<GenomicRegion::Size>(begin_distance(active_region_, read_region));
            if (diff < min_flank_padding) {
                lhs_expansion = min_flank_padding - diff;
            }
        }
        if (ends_before(active_region_, read_region)) {
            rhs_expansion = end_distance(active_region_, read_region) + min_flank_padding;
        } else {
            const auto diff = static_cast<GenomicRegion::Size>(end_distance(read_region, active_region_));
            if (diff < min_flank_padding) {
                rhs_expansion = min_flank_padding - diff;
            }
        }
        if (active_region_.begin() < lhs_expansion) {
            rhs_expansion += lhs_expansion - active_region_.begin();
            lhs_expansion = active_region_.begin();
        }
        return expand(active_region_, lhs_expansion, rhs_expansion);
    }
    return active_region_;
}

// Builder

HaplotypeGenerator::Builder& HaplotypeGenerator::Builder::set_lagging_policy(const Policies::Lagging policy) noexcept
{
    policies_.lagging = policy;
    return *this;
}

HaplotypeGenerator::Builder& HaplotypeGenerator::Builder::set_extension_policy(Policies::Extension policy) noexcept
{
    policies_.extension = policy;
    return *this;
}

HaplotypeGenerator::Builder& HaplotypeGenerator::Builder::set_target_limit(const unsigned n) noexcept
{
    policies_.haplotype_limits.target = n;
    if (policies_.haplotype_limits.holdout >= n) {
        policies_.haplotype_limits.holdout = policies_.haplotype_limits.holdout + 1;
        policies_.haplotype_limits.overflow = policies_.haplotype_limits.holdout + 1;
    }
    return *this;
}

HaplotypeGenerator::Builder& HaplotypeGenerator::Builder::set_holdout_limit(const unsigned n) noexcept
{
    policies_.haplotype_limits.holdout = std::max(n, policies_.haplotype_limits.target + 1);
    return *this;
}

HaplotypeGenerator::Builder& HaplotypeGenerator::Builder::set_overflow_limit(const unsigned n) noexcept
{
    policies_.haplotype_limits.overflow = std::max(n, policies_.haplotype_limits.holdout);
    return *this;
}

HaplotypeGenerator::Builder& HaplotypeGenerator::Builder::set_max_holdout_depth(const unsigned n) noexcept
{
    policies_.max_holdout_depth = n;
    return *this;
}

HaplotypeGenerator::Builder& HaplotypeGenerator::Builder::set_min_flank_pad(const Haplotype::MappingDomain::Size n) noexcept
{
    policies_.min_flank_pad = n;
    return *this;
}

HaplotypeGenerator::Builder& HaplotypeGenerator::Builder::set_exclusion_threshold(Haplotype::NucleotideSequence::size_type n) noexcept
{
    policies_.exclusion_threshold = n;
    return *this;
}

HaplotypeGenerator HaplotypeGenerator::Builder::build(const ReferenceGenome& reference,
                                                      const MappableFlatSet<Variant>& candidates,
                                                      const ReadMap& reads) const
{
    return HaplotypeGenerator {reference, candidates, reads, policies_};
}

} // namespace coretools
} // namespace octopus
