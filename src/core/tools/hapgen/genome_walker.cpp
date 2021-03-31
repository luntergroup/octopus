// Copyright (c) 2015-2020 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include "genome_walker.hpp"

#include <iterator>
#include <algorithm>
#include <cmath>
#include <cassert>

#include "basics/genomic_region.hpp"
#include "concepts/mappable.hpp"
#include "basics/aligned_read.hpp"
#include "core/types/variant.hpp"
#include "utils/mappable_algorithms.hpp"
#include "containers/mappable_map.hpp"

namespace octopus { namespace coretools {

GenomeWalker::GenomeWalker(Config config)
: config_ {config}
{}

namespace {

bool use_read_templates_for_lagging(const GenomeWalker::ReadTemplatePolicy policy) noexcept
{
    return policy == GenomeWalker::ReadTemplatePolicy::indicators || policy == GenomeWalker::ReadTemplatePolicy::indicators_and_extension;
}

bool use_read_templates_for_extension(const GenomeWalker::ReadTemplatePolicy policy) noexcept
{
    return policy == GenomeWalker::ReadTemplatePolicy::extension || policy == GenomeWalker::ReadTemplatePolicy::indicators_and_extension;
}

template <typename BidirIt>
BidirIt find_first_shared_helper(const ReadMap& reads, boost::optional<const TemplateMap&> read_templates,
                                 BidirIt first, BidirIt last, const Allele& allele)
{
    BidirIt result;
    if (read_templates) {
        result = find_first_shared(*read_templates, first, last, allele);
    } else {
        result = find_first_shared(reads, first, last, allele);
    }
    if (result != last) {
        while (result != first && (are_adjacent(*std::prev(result), *result) || overlaps(*std::prev(result), *result))) --result;
    }
    return result;
}

template <typename BidirIt>
bool is_sandwich_allele(BidirIt first, BidirIt allele, BidirIt last)
{
    if (allele != first && allele != last && std::next(allele) != last) {
        return has_overlapped(first, std::prev(allele), *allele) && has_overlapped(std::next(allele), last, *allele);
    } else {
        return false;
    }
}

template <typename BidirIt>
bool is_indel(const BidirIt first_itr, const BidirIt target_itr)
{
    bool result {false};
    std::find_if(std::make_reverse_iterator(target_itr), std::make_reverse_iterator(first_itr),
                 [&] (const auto& allele) {
        if (is_indel(allele)) result = true;
        return result || !is_same_region(allele, *target_itr); });
    return result;
}

template <typename BidirIt>
bool is_indel_boundary(BidirIt first, BidirIt allele, BidirIt last)
{
    if (allele != first && allele != last && std::next(allele) != last && is_indel(first, allele)) {
        auto itr = std::find_if(std::make_reverse_iterator(std::prev(allele)), std::make_reverse_iterator(first),
                                [allele] (const auto& a) { return !overlaps(a, *allele) || is_indel(a); });
        return itr != std::make_reverse_iterator(first) && overlaps(*itr, *allele) && is_indel(*itr);
    } else {
        return false;
    }
}

template <typename BidirIt>
bool is_interacting_indel(BidirIt first, BidirIt allele, BidirIt last,
                          const GenomicRegion::Size max_gap = 3)
{
    if (allele != first && allele != last && std::next(allele) != last && is_indel(first, allele)) {
        const auto interaction_region = expand_lhs(mapped_region(*allele), std::min(reference_distance(*allele), max_gap));
        auto itr = std::find_if(std::make_reverse_iterator(std::prev(allele)), std::make_reverse_iterator(first),
                                [&interaction_region] (const auto& a) { return overlaps(a, interaction_region); });
        return itr != std::make_reverse_iterator(first);
    } else {
        return false;
    }
}

template <typename BidirIt>
bool is_good_indicator_begin(BidirIt first_possible, BidirIt allele_itr, BidirIt last_possible)
{
    return !(is_sandwich_allele(first_possible, allele_itr, last_possible)
             || is_indel_boundary(first_possible, allele_itr, last_possible)
             || is_interacting_indel(first_possible, allele_itr, last_possible));
}

template <typename BidirIt>
BidirIt find_indicator_begin(BidirIt first_possible, BidirIt ideal, BidirIt last_possible)
{
    if (first_possible != ideal && ideal != last_possible) {
        auto passed_region = encompassing_region(first_possible, ideal);
        auto indicator_region = encompassing_region(ideal, last_possible);
        while (overlaps(passed_region, indicator_region)) {
            ideal = find_first_after(ideal, last_possible, passed_region);
            if (ideal == last_possible) return last_possible;
            passed_region = encompassing_region(first_possible, ideal);
            indicator_region = encompassing_region(ideal, last_possible);
        }
        while (ideal != last_possible && !is_good_indicator_begin(first_possible, ideal, last_possible)) {
            ideal = find_next_mutually_exclusive(ideal, last_possible);
        }
    }
    return ideal;
}

template <typename BidirIt>
auto get_first_included(const BidirIt first_previous_itr, const BidirIt first_included,
                        const unsigned num_indicators)
{
    auto result = first_included;
    if (num_indicators > 0) {
        assert(num_indicators <= std::distance(first_previous_itr, first_included));
        result = find_indicator_begin(first_previous_itr, std::prev(result, num_indicators), first_included);
    }
    return result;
}

template <typename BidirIt>
bool is_close(const BidirIt proposed_included, const BidirIt first_excluded)
{
    return inner_distance(*std::prev(proposed_included), *proposed_included)
           <= inner_distance(*proposed_included, *first_excluded);
}

template <typename BidirIt>
bool increases_density(const BidirIt proposed_included, const BidirIt last, const ReadMap& reads,
                       const unsigned max_density_increase)
{
    return max_count_if_shared_with_first(reads, proposed_included, last) >= max_density_increase;
}

template <typename BidirIt>
bool is_optimal_to_extend(const BidirIt first_included, const BidirIt proposed_included,
                          const BidirIt first_excluded, const BidirIt last,
                          const ReadMap& reads, const unsigned max_density_increase)
{
    if (proposed_included == last) return false;
    if (first_excluded == last) return true;
    return !increases_density(proposed_included, last, reads, max_density_increase)
           || is_close(proposed_included, first_excluded);
}

} // namespace

GenomicRegion
GenomeWalker::walk(const GenomicRegion& previous_region,
                   const ReadMap& reads,
                   const AlleleSet& alleles,
                   boost::optional<const TemplateMap&> read_templates) const
{
    using std::cbegin; using std::cend; using std::next; using std::prev; using std::min;
    using std::distance; using std::advance;
    
    if (alleles.empty()) return previous_region;
    
    auto last_allele_itr  = cend(alleles);
    auto previous_alleles = bases(overlap_range(alleles, previous_region));
    auto first_previous_itr = cbegin(previous_alleles);
    auto included_itr       = cend(previous_alleles);
    if (included_itr == last_allele_itr) {
        return shift(tail_region(rightmost_region(alleles)), 1);
    }
    if (config_.max_alleles == 0) {
        if (included_itr != last_allele_itr) {
            return *intervening_region(previous_region, *included_itr);
        } else {
            return previous_region;
        }
    }
    unsigned num_indicators {0};
    boost::optional<const TemplateMap&> indicator_read_templates {};
    if (read_templates && use_read_templates_for_lagging(config_.read_template_policy)) {
        indicator_read_templates = read_templates;
    }
    switch (config_.indicator_policy) {
        case IndicatorPolicy::includeNone: break;
        case IndicatorPolicy::includeIfSharedWithNovelRegion:
        {
            if (distance(first_previous_itr, included_itr) > 0) {
                auto it = find_first_shared_helper(reads, indicator_read_templates, first_previous_itr, included_itr, *included_itr);
                if (it != included_itr) {
                    auto expanded_leftmost = mapped_region(*it);
                    std::for_each(it, included_itr, [&] (const auto& allele) {
                        const auto ref_dist = reference_distance(allele);
                        if (ref_dist > 1) {
                            const auto max_expansion = std::min(mapped_begin(allele), 2 * ref_dist);
                            auto expanded_allele_begin = mapped_begin(allele) - max_expansion;
                            if (expanded_allele_begin < mapped_begin(expanded_leftmost)) {
                                expanded_leftmost = expand_lhs(mapped_region(allele), max_expansion);
                            }
                        }
                    });
                    for (; it != first_previous_itr && overlaps(*it, expanded_leftmost); --it);
                }
                num_indicators = static_cast<unsigned>(distance(it, included_itr));
            }
            break;
        }
        case IndicatorPolicy::includeIfLinkableToNovelRegion:
        {
            if (distance(first_previous_itr, included_itr) > 0) {
                auto it = included_itr;
                while (true) {
                    const auto it2 = find_first_shared_helper(reads, indicator_read_templates, first_previous_itr, it, *it);
                    if (it2 == it) {
                        it = it2;
                        break;
                    } else {
                        it = it2;
                    }
                }
                num_indicators = static_cast<unsigned>(distance(it, included_itr));
            }
            break;
        }
        case IndicatorPolicy::includeAll:
            num_indicators = static_cast<unsigned>(distance(first_previous_itr, included_itr));
            break;
    }
    auto first_included_itr = get_first_included(first_previous_itr, included_itr, num_indicators);
    auto num_remaining_alleles = static_cast<unsigned>(distance(included_itr, last_allele_itr));
    unsigned num_excluded_alleles {0};
    auto num_included = config_.max_alleles;
    if (config_.extension_policy == ExtensionPolicy::includeIfWithinReadLengthOfFirstIncluded) {
        auto max_alleles_within_read_length = static_cast<unsigned>(max_count_if_shared_with_first(reads, first_included_itr, last_allele_itr));
        num_included = min({num_included, num_remaining_alleles, max_alleles_within_read_length + 1});
        num_excluded_alleles = max_alleles_within_read_length - num_included;
    } else {
        num_included = min(num_included, num_remaining_alleles);
    }
    boost::optional<const TemplateMap&> extension_read_templates {};
    if (read_templates && use_read_templates_for_extension(config_.read_template_policy)) {
        extension_read_templates = read_templates;
    }
    assert(num_included > 0);
    auto first_excluded_itr = next(included_itr, num_included);
    while (--num_included > 0 && is_optimal_to_extend(first_included_itr, next(included_itr), first_excluded_itr,
                                                      last_allele_itr, reads, num_included + num_excluded_alleles)) {
        if (!can_extend(*included_itr, *next(included_itr), reads, extension_read_templates)) {
            break;
        }
        ++included_itr;
    }
    const auto rightmost = rightmost_mappable(first_included_itr, next(included_itr));
    const auto first_exclusive = find_next_mutually_exclusive(rightmost, last_allele_itr);
    return encompassing_region(first_included_itr, first_exclusive);
}

bool
GenomeWalker::can_extend(const Allele& active, const Allele& novel,
                         const ReadMap& reads, boost::optional<const TemplateMap&> read_templates) const
{
    if (config_.max_extension && inner_distance(active, novel) > *config_.max_extension) {
        return false;
    }
    if (config_.extension_policy == ExtensionPolicy::includeIfAllSamplesSharedWithFrontier) {
        if (read_templates) {
            return all_shared(*read_templates, active, novel);
        } else {
            return all_shared(reads, active, novel);
        }
    } else if (config_.extension_policy == ExtensionPolicy::includeIfAnySampleSharedWithFrontier) {
        if (read_templates) {
            return has_shared(*read_templates, active, novel);
        } else {
            return has_shared(reads, active, novel);
        }
    }
    return config_.extension_policy == ExtensionPolicy::noLimit;
}

} // namespace coretools
} // namespace octopus
