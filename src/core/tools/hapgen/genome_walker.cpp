// Copyright (c) 2017 Daniel Cooke
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

GenomeWalker::GenomeWalker(unsigned max_included,
                           IndicatorPolicy indicator_policy,
                           ExtensionPolicy extension_policy)
: max_included_ {max_included}
, indicator_policy_ {indicator_policy}
, extension_policy_ {extension_policy}
{}

GenomicRegion GenomeWalker::walk(const GenomicRegion::ContigName& contig,
                                 const ReadMap& reads,
                                 const AlleleSet& alleles) const
{
    return walk(GenomicRegion {contig, 0, 0}, reads, alleles);
}

namespace {

template <typename BidirIt>
bool is_sandwich_allele(BidirIt first, BidirIt allele, BidirIt last)
{
    if (allele != first && allele != last && std::next(allele) != last) {
        return overlaps(*std::prev(allele), *allele) && overlaps(*allele, *std::next(allele));
    } else {
        return false;
    }
}

template <typename BidirIt>
bool is_indel_boundary(BidirIt first, BidirIt allele, BidirIt last)
{
    if (allele != first && allele != last && std::next(allele) != last && is_indel(*allele)) {
        auto itr = std::find_if(std::make_reverse_iterator(std::prev(allele)), std::make_reverse_iterator(first),
                                [&] (const auto& a) { return !overlaps(a, *allele) || is_indel(a); });
        return itr != std::make_reverse_iterator(first) && overlaps(*itr, *allele) && is_indel(*itr);
    } else {
        return false;
    }
}

template <typename BidirIt>
bool is_good_indicator_begin(BidirIt first_possible, BidirIt allele_itr, BidirIt last_possible)
{
    return !(is_sandwich_allele(first_possible, allele_itr, last_possible)
             || is_indel_boundary(first_possible, allele_itr, last_possible));
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

}

GenomicRegion GenomeWalker::walk(const GenomicRegion& previous_region,
                                 const ReadMap& reads,
                                 const AlleleSet& alleles) const
{
    using std::cbegin; using std::cend; using std::next; using std::prev; using std::min;
    using std::distance; using std::advance;
    
    auto last_allele_itr  = cend(alleles);
    auto previous_alleles = bases(overlap_range(alleles, previous_region));
    auto first_previous_itr = cbegin(previous_alleles);
    auto included_itr       = cend(previous_alleles);
    if (included_itr == last_allele_itr) {
        return shift(tail_region(previous_region), 2);
    }
    if (max_included_ == 0) {
        if (included_itr != last_allele_itr) {
            return *intervening_region(previous_region, *included_itr);
        } else {
            return previous_region;
        }
    }
    unsigned num_indicators {0};
    switch (indicator_policy_) {
        case IndicatorPolicy::includeNone: break;
        case IndicatorPolicy::includeIfSharedWithNovelRegion:
        {
            if (distance(first_previous_itr, included_itr) > 0) {
                auto it = find_first_shared(reads, first_previous_itr, included_itr, *included_itr);
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
                    const auto it2 = find_first_shared(reads, first_previous_itr, it, *it);
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
    auto num_included = max_included_;
    if (extension_policy_ == ExtensionPolicy::includeIfWithinReadLengthOfFirstIncluded) {
        auto max_alleles_within_read_length = static_cast<unsigned>(max_count_if_shared_with_first(reads, first_included_itr, last_allele_itr));
        num_included = min({num_included, num_remaining_alleles, max_alleles_within_read_length + 1});
        num_excluded_alleles = max_alleles_within_read_length - num_included;
    } else {
        num_included = min(num_included, num_remaining_alleles);
    }
    assert(num_included > 0);
    auto first_excluded_itr = next(included_itr, num_included);
    while (--num_included > 0 && is_optimal_to_extend(first_included_itr, next(included_itr), first_excluded_itr,
                                                      last_allele_itr, reads, num_included + num_excluded_alleles)) {
        if (extension_policy_ == ExtensionPolicy::includeIfAllSamplesSharedWithFrontier
            && !all_shared(reads, *included_itr, *next(included_itr))) {
            break;
        } else if (extension_policy_ == ExtensionPolicy::includeIfAnySampleSharedWithFrontier
            && !has_shared(reads, *included_itr, *next(included_itr))) {
            break;
        }
        ++included_itr;
    }
    const auto rightmost = rightmost_mappable(first_included_itr, next(included_itr));
    const auto first_exclusive = find_next_mutually_exclusive(rightmost, last_allele_itr);
    return encompassing_region(first_included_itr, first_exclusive);
}

} // namespace coretools
} // namespace octopus
