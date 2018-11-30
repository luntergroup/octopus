// Copyright (c) 2015-2018 Daniel Cooke
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

auto get_contig(const MappableFlatSet<Variant>& candidates)
{
    if (!candidates.empty()) {
        return contig_name(candidates.front());
    } else {
        throw std::runtime_error {"HaplotypeGenerator: not supplied with any candidates"};
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

auto make_lagged_walker(const HaplotypeGenerator::Policies& policies)
{
    return GenomeWalker {
        max_included(policies.haplotype_limits.target),
        get_walker_policy(policies.lagging),
        get_walker_policy(policies.extension)
    };
}

bool all_empty(const ReadMap& reads) noexcept
{
    return std::all_of(std::cbegin(reads), std::cend(reads), [] (const auto& p) noexcept { return p.second.empty(); });
}

} // namespace

// public members

HaplotypeGenerator::HaplotypeGenerator(const ReferenceGenome& reference,
                                       const MappableFlatSet<Variant>& candidates,
                                       const ReadMap& reads,
                                       boost::optional<const ReadPipe::Report&> reads_report,
                                       Policies policies,
                                       DenseVariationDetector dense_variation_detector)
: policies_ {std::move(policies)}
, tree_ {get_contig(candidates), reference}
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
    assert(!candidates.empty());
    assert(!alleles_.empty());
    if (policies.lagging != Policies::Lagging::none) {
        lagged_walker_ = make_lagged_walker(policies);
    }
    if (!all_empty(reads_)) {
        const auto dense_regions = dense_variation_detector.detect(candidates, reads, reads_report);
        for (const auto& dense : dense_regions) {
            if (dense.action == DenseVariationDetector::DenseRegion::RecommendedAction::skip) {
                if (debug_log_) {
                    stream(*debug_log_) << "Erasing " << count_contained(alleles_, dense.region)
                                        << " alleles in dense region " << dense.region;
                }
                alleles_.erase_contained(dense.region);
            } else if (is_lagging_enabled()) {
                lagging_exclusion_zones_.insert(dense.region);
            }
        }
        if (!lagging_exclusion_zones_.empty() && debug_log_) {
            auto log = stream(*debug_log_);
            log << "Found lagging exclusion zones: ";
            for (const auto& zone : lagging_exclusion_zones_) log << zone << " ";
        }
        if (alleles_.empty()) {
            alleles_.insert(candidates.back().ref_allele());
        }
    }
    rightmost_allele_ = alleles_.rightmost();
    active_region_ = head_region(alleles_.leftmost());
    if (active_region_.begin() != 0) {
        active_region_ = shift(active_region_, -1);
    }
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
    if (alleles_.empty()) {
        return std::make_tuple(std::vector<Haplotype> {}, boost::none, boost::none);
    }
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
    GenomicRegion passed_region {active_region_.contig_name(), 0, active_region_.begin()};
    alleles_.erase_contained(passed_region);
}

void HaplotypeGenerator::jump(GenomicRegion region)
{
    clear_progress();
    next_active_region_ = std::move(region);
    remove_passed_alleles();
}

bool HaplotypeGenerator::removal_has_impact() const
{
    if (in_holdout_mode()) return true;
    if (!is_lagging_enabled(active_region_) || contains(active_region_, rightmost_allele_)) return false;
    const auto max_lagged_region = lagged_walker_->walk(active_region_, reads_, alleles_);
    return overlaps(max_lagged_region, active_region_);
}

unsigned HaplotypeGenerator::max_removal_impact() const
{
    if (in_holdout_mode()) return tree_.num_haplotypes();
    if (!is_lagging_enabled(active_region_) || contains(active_region_, rightmost_allele_)) return 0;
    const auto max_lagged_region = lagged_walker_->walk(active_region_, reads_, alleles_);
    if (!overlaps(max_lagged_region, active_region_)) return 0;
    const auto novel_region = right_overhang_region(max_lagged_region, active_region_);
    const auto num_novel_alleles = count_overlapped(alleles_, novel_region);
    if (num_novel_alleles == 0) return 0;
    unsigned num_leftover_haplotypes {0};
    static auto max_exponent = static_cast<unsigned>(std::log2(std::numeric_limits<unsigned>::max()));
    if (num_novel_alleles / 2 < max_exponent) {
        const auto max_new_haplotypes = std::max(static_cast<unsigned>(std::exp2(num_novel_alleles / 2)), 1u);
        num_leftover_haplotypes = policies_.haplotype_limits.target / max_new_haplotypes;
    }
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

bool HaplotypeGenerator::is_lagging_enabled(const GenomicRegion& region) const
{
    return is_lagging_enabled() && !has_overlapped(lagging_exclusion_zones_, region);
}

bool HaplotypeGenerator::is_active_region_lagged() const
{
    if (in_holdout_mode()) return true;
    if (!is_lagging_enabled(active_region_)) return false;
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
        if (is_lagging_enabled(active_region_) || in_holdout_mode()) {
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

bool has_duplicate_insertion(const std::vector<GenomicRegion>& indicator_blocks, const std::vector<GenomicRegion>& novel_blocks)
{
    return is_empty(novel_blocks.front()) && are_adjacent(indicator_blocks.back(), novel_blocks.front());
}

auto max_ref_distance(const GenomicRegion& region, const MappableFlatSet<Allele>& alleles)
{
    const auto contained = contained_range(alleles, region);
    auto itr = std::max_element(std::cbegin(contained), std::cend(contained),
                                [] (const auto& lhs, const auto& rhs) {
                                    return reference_distance(lhs) < reference_distance(rhs);
                                });
    return reference_distance(*itr);
}

auto max_ref_distance(const GenomicRegion& region, const HaplotypeTree& tree)
{
    HaplotypeTree::HaplotypeLength min_length, max_length;
    std::tie(min_length, max_length) = minmax_haplotype_lengths(tree, region);
    if (min_length == max_length) {
        return HaplotypeTree::HaplotypeLength {0};
    } else {
        if (min_length <= size(region) && size(region) <= max_length) {
            return std::max(max_length - size(region), size(region) - min_length);
        } else if (min_length <= size(region)) {
            return size(region) - min_length;
        } else if (size(region) <= max_length) {
            return max_length - size(region);
        } else {
            return HaplotypeTree::HaplotypeLength {0};
        }
    }
}

GenomicRegion::Size max_ref_distance(const GenomicRegion& region, const MappableFlatSet<Allele>& alleles, const HaplotypeTree& tree)
{
    if (max_ref_distance(region, alleles) > 0) {
        return max_ref_distance(region, tree);
    } else {
        return 0;
    }
}

std::vector<unsigned> get_max_ref_distances(const std::vector<GenomicRegion>& regions, const MappableFlatSet<Allele>& alleles,
                                            const HaplotypeTree& tree)
{
    std::vector<unsigned> result(regions.size());
    std::transform(std::cbegin(regions), std::cend(regions), std::begin(result),
                   [&] (const auto& region) { return max_ref_distance(region, alleles, tree); });
    return result;
}

std::vector<GenomicRegion> get_joined_blocks(const std::vector<GenomicRegion>& blocks,
                                             const MappableFlatSet<Allele>& alleles,
                                             const HaplotypeTree& tree)
{
    std::vector<GenomicRegion> expanded_blocks {};
    expanded_blocks.reserve(blocks.size());
    for (const auto& block : blocks) {
        auto ref_dist = max_ref_distance(block, alleles, tree);
        if (ref_dist > 0) expanded_blocks.push_back(expand(block, ref_dist));
    }
    if (expanded_blocks.empty()) return blocks;
    std::sort(std::begin(expanded_blocks), std::end(expanded_blocks));
    std::vector<GenomicRegion> interacting_expanded_blocks {};
    interacting_expanded_blocks.reserve(expanded_blocks.size());
    std::copy_if(std::begin(expanded_blocks), std::end(expanded_blocks), std::back_inserter(interacting_expanded_blocks),
                 [&] (const auto& block) { return count_overlapped(expanded_blocks, block) > 1; });
    if (interacting_expanded_blocks.empty()) return blocks;
    const auto join_regions = extract_covered_regions(interacting_expanded_blocks);
    std::vector<GenomicRegion> result {};
    result.reserve(blocks.size());
    auto last_block = std::cbegin(blocks);
    for (auto&& join_region : join_regions) {
        auto joined = bases(contained_range(last_block, std::cend(blocks), join_region));
        result.insert(std::cend(result), last_block, std::begin(joined));
        result.push_back(closed_region(joined.front(), joined.back()));
        last_block = std::end(joined);
    }
    result.insert(std::cend(result), last_block, std::cend(blocks));
    return result;
}

bool requires_exclusion_zone(const GenomicRegion& region, const MappableFlatSet<Allele>& alleles,
                             const HaplotypeTree& tree, const Haplotype::NucleotideSequence::size_type ref_distance_join_threshold)
{
    return max_ref_distance(region, alleles) > ref_distance_join_threshold
           && max_ref_distance(region, tree) > ref_distance_join_threshold;
}

boost::optional<GenomicRegion>
get_exclusion_zone(const GenomicRegion& region, const MappableFlatSet<Allele>& alleles,
                   const HaplotypeTree& tree, const Haplotype::NucleotideSequence::size_type ref_distance_join_threshold)
{
    if (requires_exclusion_zone(region, alleles, tree, ref_distance_join_threshold)) {
        return expand(region, 2 * ref_distance_join_threshold);
    } else {
        return boost::none;
    }
}

std::vector<GenomicRegion>
get_exclusion_zones(const std::vector<GenomicRegion>& regions, const MappableFlatSet<Allele>& alleles,
                    const HaplotypeTree& tree, const Haplotype::NucleotideSequence::size_type ref_distance_join_threshold)
{
    std::vector<GenomicRegion> result {};
    result.reserve(regions.size());
    for (const auto& region : regions) {
        auto exclusion_zone = get_exclusion_zone(region, alleles, tree, ref_distance_join_threshold);
        if (exclusion_zone) result.push_back(std::move(*exclusion_zone));
    }
    if (result.size() > 1) {
        result = extract_mutually_exclusive_regions(result);
    }
    return result;
}

auto merge(const std::vector<GenomicRegion>& indicator_blocks, const std::vector<GenomicRegion>& exclusion_zones)
{
    std::vector<GenomicRegion> result {};
    result.reserve(indicator_blocks.size());
    auto last_block = std::cbegin(indicator_blocks);
    for (auto&& exclusion_zone : exclusion_zones) {
        auto interacting = bases(overlap_range(last_block, std::cend(indicator_blocks), exclusion_zone, BidirectionallySortedTag {}));
        result.insert(std::cend(result), last_block, std::begin(interacting));
        if (empty(interacting)) {
            result.push_back(std::move(exclusion_zone));
        } else {
            auto indicator_block = closed_region(interacting.front(), interacting.back());
            result.push_back(encompassing_region(indicator_block, exclusion_zone));
        }
        last_block = std::end(interacting);
    }
    result.insert(std::cend(result), last_block, std::cend(indicator_blocks));
    return result;
}

auto extract_indicator_alleles(std::vector<GenomicRegion> indicator_blocks, const MappableFlatSet<Allele>& alleles,
                               const HaplotypeTree& tree,
                               const boost::optional<Haplotype::NucleotideSequence::size_type> max_indicator_join_distance)
{
    if (indicator_blocks.size() < 2) return indicator_blocks;
    indicator_blocks = get_joined_blocks(std::move(indicator_blocks), alleles, tree);
    if (!max_indicator_join_distance) return indicator_blocks;
    auto exclusion_zones = get_exclusion_zones(indicator_blocks, alleles, tree, *max_indicator_join_distance);
    if (exclusion_zones.empty()) {
        return indicator_blocks;
    } else {
        return merge(std::move(indicator_blocks), std::move(exclusion_zones));
    }
}

double get_required_novel_fraction(HaplotypeGenerator::Policies::Lagging policy) noexcept
{
    using HGP = HaplotypeGenerator::Policies::Lagging;
    switch (policy) {
        case HGP::none:
        case HGP::conservative:
            return 1.5;
        case HGP::moderate:
            return 2;
        case HGP::normal:
            return 4;
        case HGP::aggressive:
            return 6;
        default:
            return 1; // prevents compiler warning
    }
}

auto to_log2s(const std::vector<unsigned>& values)
{
    std::vector<double> result(values.size());
    std::transform(std::cbegin(values), std::cend(values), std::begin(result),
                   [] (auto value) { return std::log2(value); });
    return result;
}

auto extract_overlap_counts(const std::vector<GenomicRegion>& regions, const MappableFlatSet<Allele>& alleles)
{
    std::vector<unsigned> result(regions.size());
    std::transform(std::cbegin(regions), std::cend(regions), std::begin(result),
                   [&] (const auto& region) { return count_overlapped(alleles, region); });
    return result;
}

auto get_partial_overlap_factors(const std::vector<GenomicRegion>& regions, const MappableFlatSet<Allele>& alleles)
{
    auto factors = to_log2s(extract_overlap_counts(regions, alleles));
    std::partial_sum(std::cbegin(factors), std::cend(factors), std::begin(factors));
    std::vector<std::size_t> result(factors.size());
    std::transform(std::cbegin(factors), std::cend(factors), std::begin(result),
                   [] (auto f) { return static_cast<std::size_t>(std::round(f)); });
    return result;
}

unsigned lagging_level(const HaplotypeGenerator::Policies::Lagging& policy) noexcept
{
    using HGP = HaplotypeGenerator::Policies::Lagging;
    switch (policy) {
        case HGP::none: return 0;
        case HGP::conservative: return 1;
        case HGP::moderate: return 2;
        case HGP::normal: return 3;
        case HGP::aggressive: return 4;
        default:
            return 0; // prevents compiler warning
    }
}

std::size_t get_target_tree_size(const HaplotypeTree curr_tree,
                                 const std::vector<GenomicRegion>& indicator_blocks,
                                 const std::vector<GenomicRegion>& novel_blocks,
                                 const MappableFlatSet<Allele>& alleles,
                                 const HaplotypeGenerator::Policies policies)
{
    const auto curr_tree_size = curr_tree.num_haplotypes();
    if (curr_tree_size == 0) return 0;
    if (novel_blocks.empty() || indicator_blocks.empty()) return curr_tree_size;
    const auto novel_factors = get_partial_overlap_factors(novel_blocks, alleles);
    auto target_tree_size = policies.haplotype_limits.target;
    const auto effective_log_space = static_cast<std::size_t>(std::log2(std::max(target_tree_size / curr_tree_size, std::size_t {1})));
    if (effective_log_space >= novel_factors.back()) {
        return curr_tree_size;
    } else {
        const auto novel_fraction = get_required_novel_fraction(policies.lagging);
        const auto num_novels_required = static_cast<std::size_t>(novel_factors.size() / novel_fraction);
        auto required_novel_factor = novel_factors[num_novels_required > 0 ? num_novels_required - 1 : 0];
        if (novel_blocks.size() == 1 && lagging_level(policies.lagging) > 1) {
            target_tree_size = policies.haplotype_limits.holdout;
        }
        const auto tree_reduction_denom = static_cast<std::size_t>(std::pow(2, required_novel_factor));
        return std::min(target_tree_size / tree_reduction_denom, curr_tree_size);
    }
}

auto expand_rhs_by_max_ref_dist(const GenomicRegion& region, const MappableFlatSet<Allele>& alleles, const HaplotypeTree tree)
{
    return expand_rhs(region,  max_ref_distance(region, alleles, tree));
}

auto find_rightmost_expanded(const std::vector<GenomicRegion>& blocks, const MappableFlatSet<Allele>& alleles,
                             const HaplotypeTree tree)
{
    std::vector<GenomicRegion> expanded_blocks(blocks.size());
    std::transform(std::cbegin(blocks), std::cend(blocks), std::begin(expanded_blocks),
                   [&] (const auto& block) { return expand_rhs_by_max_ref_dist(block, alleles, tree); });
    const auto itr = rightmost_mappable(expanded_blocks);
    const auto offset =  std::distance(std::cbegin(expanded_blocks), itr);
    return std::make_pair(std::next(std::cbegin(blocks), offset), *itr);
}

auto expand_lhs_by_max_ref_dist(const GenomicRegion& region, const MappableFlatSet<Allele>& alleles, const HaplotypeTree tree)
{
    const auto max_ref_dist = max_ref_distance(region, alleles, tree);
    return expand_lhs(region, std::min(max_ref_dist, static_cast<std::remove_const_t<decltype(max_ref_dist)>>(region.begin())));
}

auto get_leftmost_expanded(const std::vector<GenomicRegion>& blocks, const MappableFlatSet<Allele>& alleles,
                           const HaplotypeTree tree)
{
    assert(!blocks.empty());
    auto result = expand_lhs_by_max_ref_dist(blocks.front(), alleles, tree);
    std::for_each(std::next(std::cbegin(blocks)), std::cend(blocks), [&] (const auto& block) {
        auto expanded_block = expand_lhs_by_max_ref_dist(block, alleles, tree);
        if (begins_before(expanded_block, result)) {
            result = std::move(expanded_block);
        }
    });
    return result;
}

std::vector<GenomicRegion>
remove_interacting_indicator_tail(std::vector<GenomicRegion>& indicator_blocks,
                               const std::vector<GenomicRegion>& novel_blocks,
                               const MappableFlatSet<Allele>& alleles, const HaplotypeTree tree)
{
    std::vector<GenomicRegion> result {};
    if (!indicator_blocks.empty()) {
        std::vector<GenomicRegion>::const_iterator rightmost_indicator_itr;
        GenomicRegion expanded_indicator_tail;
        std::tie(rightmost_indicator_itr, expanded_indicator_tail) = find_rightmost_expanded(indicator_blocks, alleles, tree);
        const auto leftmost_novel = get_leftmost_expanded(novel_blocks, alleles, tree);
        if (overlaps(expanded_indicator_tail, leftmost_novel)) {
            result.assign(rightmost_indicator_itr, std::cend(indicator_blocks));
            indicator_blocks.erase(rightmost_indicator_itr, std::cend(indicator_blocks));
        }
    }
    return result;
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

void remove_blocks_after(std::vector<GenomicRegion>& blocks, const GenomicRegion& region)
{
    auto itr = std::find_if_not(std::crbegin(blocks), std::crend(blocks),
                                [&] (const auto& block) { return is_after(block, region); });
    blocks.erase(itr.base(), std::cend(blocks));
}

unsigned get_num_ideal_new_novel_blocks(const std::vector<GenomicRegion>& novel_blocks,
                                        const GenomicRegion& indicator_region,
                                        const MappableFlatSet<Allele>& alleles)
{
    const auto num_indicators = count_contained(alleles, indicator_region);
    if (num_indicators > 20) {
        if (novel_blocks.size() > 10) {
            return 4;
        } else if (novel_blocks.size() > 5) {
            return 3;
        } else if (novel_blocks.size() >= 3) {
            return 2;
        }
    }
    return 1;
}

template <typename Range>
unsigned extend_novel(HaplotypeTree& tree, const std::vector<GenomicRegion>& novel_blocks,
                      const Range& novel_alleles, const unsigned ideal_num_added_novel_blocks,
                      const HaplotypeGenerator::Policies::HaplotypeLimits& limits)
{
    unsigned num_novel_blocks_added {0};
    for (const auto& block : novel_blocks) {
        const auto interacting_alleles = contained_range(novel_alleles, block);
        const auto last_added = extend_tree_until(interacting_alleles, tree, limits.overflow);
        if (last_added != std::cend(interacting_alleles)) {
            if (num_novel_blocks_added > 0) {
                tree.clear(block);
            }
            break;
        }
        ++num_novel_blocks_added;
        if (num_novel_blocks_added >= ideal_num_added_novel_blocks || tree.num_haplotypes() >= limits.holdout) {
            if (tree.num_haplotypes() > limits.target) {
                if (num_novel_blocks_added > ideal_num_added_novel_blocks
                    && !(num_novel_blocks_added == novel_blocks.size() && tree.num_haplotypes() < limits.holdout)) {
                    tree.clear(block);
                    --num_novel_blocks_added;
                }
                break;
            } else if (tree.num_haplotypes() == limits.target) {
                break;
            }
        }
    }
    return num_novel_blocks_added;
}

} // namespace

void HaplotypeGenerator::update_lagged_next_active_region() const
{
    if (contains(active_region_, rightmost_allele_)) {
        next_active_region_ = shift(tail_region(rightmost_allele_), 2); // Nothing more to do
        return;
    }
    const auto max_lagged_region = find_max_lagged_region();
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
        auto novel_blocks = extract_mutually_exclusive_regions(novel_alleles);
        // contained_range is required for indicators as previous active alleles must be contained
        // within the reported active region.
        const auto indicator_region = *overlapped_region(active_region_, max_lagged_region);
        std::vector<GenomicRegion> protected_indicator_blocks {};
        std::size_t target_tree_size {};
        if (has_contained(alleles_, indicator_region)) {
            const auto indicator_alleles = contained_range(alleles_, indicator_region);
            auto indicator_blocks = extract_mutually_exclusive_regions(indicator_alleles);
            // Although HaplotypeTree can handle duplicate extension, if the tree is not in a complete
            // state (i.e. it has been pruned) then adding duplicates will corrupt the tree.
            if (has_duplicate_insertion(indicator_blocks, novel_blocks)) {
                pop_front(novel_blocks);
            }
            if (!in_holdout_mode()) {
                const auto novel_overlap_region = encompassing_region(novel_alleles);
                while (!indicator_blocks.empty() && overlaps(indicator_blocks.back(), novel_overlap_region)) {
                    // This can only happen if holdouts have just been reintroduced, in which case we may need to
                    // 'rewind' some indicators
                    indicator_blocks.pop_back();
                }
                auto expanded_indicator_blocks = extract_indicator_alleles(indicator_blocks, alleles_, tree_,
                                                                           policies_.max_indicator_join_distance);
                target_tree_size = get_target_tree_size(test_tree, expanded_indicator_blocks, novel_blocks,
                                                        alleles_, policies_);
                protected_indicator_blocks = remove_interacting_indicator_tail(expanded_indicator_blocks, novel_blocks,
                                                                               alleles_, tree_);
                prune_indicators(test_tree, expanded_indicator_blocks, target_tree_size);
            }
        }
        if (in_holdout_mode() && ends_before(top_holdout_region(), novel_region)) {
            remove_blocks_after(novel_blocks, top_holdout_region());
        }
        const auto ideal_num_new_novel_blocks = get_num_ideal_new_novel_blocks(novel_blocks, indicator_region, alleles_);
        auto num_novel_blocks_added = extend_novel(test_tree, novel_blocks, novel_alleles, ideal_num_new_novel_blocks,
                                                   policies_.haplotype_limits);
        if (num_novel_blocks_added == 0 && protected_indicator_blocks.size() > 1) {
            auto last_block = protected_indicator_blocks.back();
            protected_indicator_blocks.pop_back();
            prune_indicators(test_tree, protected_indicator_blocks, target_tree_size);
            protected_indicator_blocks.assign({last_block});
            num_novel_blocks_added = extend_novel(test_tree, novel_blocks, novel_alleles, 1, policies_.haplotype_limits);
        }
        if (num_novel_blocks_added > 0) {
            next_active_region_ = test_tree.encompassing_region();
        } else {
            assert(!novel_blocks.empty());
            if (test_tree.is_empty()) {
                next_active_region_ = novel_blocks.front();
            } else {
                next_active_region_ = closed_region(test_tree.encompassing_region(), novel_blocks.front());
            }
        }
    }
}

void HaplotypeGenerator::remove_passed_alleles()
{
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
    if (!in_holdout_mode()) remove_passed_alleles();
    auto novel_active_region = *next_active_region_;
    if (!tree_.is_empty()) {
        novel_active_region = right_overhang_region(*next_active_region_, active_region_);
    }
    auto novel_active_alleles = overlap_range(alleles_, novel_active_region);
    assert(!empty(novel_active_alleles));
    if (overlaps(novel_active_alleles.front(), active_region_)) {
        // Then there are duplicate insertions at the front of the novel region that must be removed
        while (overlaps(novel_active_alleles.front(), active_region_)) {
            novel_active_alleles.advance_begin(1);
            assert(!empty(novel_active_alleles));
        }
        novel_active_region = encompassing_region(novel_active_alleles);
    }
    auto last_added_novel_itr = extend_tree_until(novel_active_alleles, tree_, policies_.haplotype_limits.holdout);
    if (last_added_novel_itr != std::cend(novel_active_alleles)) {
        reset_next_active_region();
        const auto active_region_before_holdout = active_region_;
        const auto novel_region_before_holdout = novel_active_region;
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
            if (begins_before(novel_active_region, novel_region_before_holdout)) {
                if (!is_before(novel_active_region, novel_region_before_holdout)) {
                    novel_active_region = closed_region(novel_region_before_holdout, novel_active_region);
                }
            }
            novel_active_alleles = overlap_range(alleles_, novel_active_region);
            assert(!empty(novel_active_alleles));
            while (begins_before(novel_active_alleles.front(), active_region_before_holdout)
                   && ends_before(novel_active_alleles.front(), novel_active_region)) {
                novel_active_alleles.advance_begin(1);
                assert(!empty(novel_active_alleles));
            }
            last_added_novel_itr = extend_tree_until(novel_active_alleles, tree_, policies_.haplotype_limits.holdout);
            if (overlaps(active_region_, top_holdout_region())) {
                next_holdout_region = *overlapped_region(active_region_, top_holdout_region());
            } else {
                next_holdout_region = novel_active_region;
            }
        }
        if (last_added_novel_itr != std::cend(novel_active_alleles)) {
            last_added_novel_itr = extend_tree_until(last_added_novel_itr, std::cend(novel_active_alleles), tree_,
                                                     policies_.haplotype_limits.overflow);
            if (last_added_novel_itr != std::cend(novel_active_alleles)) {
                if (in_holdout_mode()) {
                    active_region_ = encompassing_region(active_region_, tree_.encompassing_region());
                    active_region_ = encompassing_region(active_region_, *holdout_region_);
                } else {
                    active_region_ = tree_.encompassing_region();
                }
                throw HaplotypeOverflow {active_region_, tree_.num_haplotypes()};
            } else {
                active_region_ = tree_.encompassing_region();
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

auto extract_unique_regions(const MappableFlatSet<Allele>& alleles)
{
    auto result = extract_regions(alleles);
    result.erase(std::unique(std::begin(result), std::end(result)), std::end(result));
    return result;
}

struct ContainmentCount
{
    GenomicRegion region;
    std::size_t count;
};

auto get_containment_counts(const MappableFlatSet<Allele>& alleles)
{
    const auto allele_regions = extract_unique_regions(alleles);
    std::vector<ContainmentCount> result {};
    result.reserve(result.size());
    std::transform(std::cbegin(allele_regions), std::cend(allele_regions), std::back_inserter(result),
                   [&alleles] (const auto& region) -> ContainmentCount {
                       return {region, count_fully_contained(alleles, region)};
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
    const auto active_alleles = copy_overlapped(alleles_, region);
    assert(!active_alleles.empty());
    auto containment_counts = get_containment_counts(active_alleles);
    auto last_viable_holdout = std::upper_bound(std::begin(containment_counts), std::end(containment_counts),
                                                3u, [] (auto lhs, const auto& rhs) { return lhs > rhs.count; });
    if (last_viable_holdout != std::end(containment_counts)) {
        last_viable_holdout = std::remove_if(std::begin(containment_counts), last_viable_holdout,
                                             [] (const auto& c) { return size(c.region) < 2; });
        const auto num_viable_holdouts = std::distance(std::begin(containment_counts), last_viable_holdout);
        if (num_viable_holdouts == 0) return boost::none;
        auto holdout_itr = std::begin(containment_counts);
        if (num_viable_holdouts > 1) {
            const auto good_count = std::log2(policies_.haplotype_limits.holdout);
            auto itr = std::find_if(std::next(holdout_itr), last_viable_holdout,
                                    [=] (const auto& c) {
                                        return c.count > good_count && is_before(c.region, holdout_itr->region);
                                    });
            if (itr != last_viable_holdout) holdout_itr = itr;
        }
        const auto holdouts = get_holdout_range(active_alleles, holdout_itr->region);
        assert(!holdouts.empty());
        return HoldoutSet {std::cbegin(holdouts), std::cend(holdouts), holdout_itr->region};
    } else {
        return boost::none;
    }
}

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
auto sum_indel_sizes(const Range& alleles) noexcept
{
    return std::accumulate(std::cbegin(alleles), std::cend(alleles), std::size_t {0},
                           [] (const auto curr, const Allele& allele) noexcept {
                               return curr + reference_distance(allele);
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
    return expand(active_region_, 1);
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

HaplotypeGenerator::Builder& HaplotypeGenerator::Builder::set_max_indicator_join_distance(Haplotype::NucleotideSequence::size_type n) noexcept
{
    policies_.max_indicator_join_distance = n;
    return *this;
}

HaplotypeGenerator::Builder& HaplotypeGenerator::Builder::set_max_expected_log_allele_count_per_base(double v) noexcept
{
    policies_.max_expected_log_allele_count_per_base = v;
    return *this;
}

HaplotypeGenerator::Builder& HaplotypeGenerator::Builder::set_dense_variation_detector(DenseVariationDetector detector) noexcept
{
    dense_variation_detector_ = std::move(detector);
    return *this;
}

HaplotypeGenerator HaplotypeGenerator::Builder::build(const ReferenceGenome& reference,
                                                      const MappableFlatSet<Variant>& candidates,
                                                      const ReadMap& reads,
                                                      boost::optional<const ReadPipe::Report&> reads_report) const
{
    return HaplotypeGenerator {reference, candidates, reads, std::move(reads_report), policies_, dense_variation_detector_};
}

} // namespace coretools
} // namespace octopus
