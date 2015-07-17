//
//  search_regions.h
//  Octopus
//
//  Created by Daniel Cooke on 27/04/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_search_regions_h
#define Octopus_search_regions_h

#include <iterator>  // std::distance, std::next, std::prev, std::advance, std::cbegin, std::cend
#include <algorithm> // std::min
#include <cmath>     // std::abs

#include "genomic_region.h"
#include "region_algorithms.h"
#include "map_utils.h"

#include <iostream> // TEST

namespace Octopus
{

namespace detail
{
    template <typename BidirectionalIterator, typename SampleReadMap>
    bool is_optimal_to_extend(BidirectionalIterator first_included, BidirectionalIterator proposed_included,
                              BidirectionalIterator first_excluded, BidirectionalIterator last,
                              const SampleReadMap& reads, unsigned max_density_increase)
    {
        if (proposed_included == last) return false;
        if (first_excluded == last) return true;
        
        bool increases_density {max_count_if_shared_with_first(reads, proposed_included, last)
                                        >= max_density_increase};
        
        return (increases_density) ? inner_distance(*std::prev(proposed_included), *proposed_included)
                                        <= inner_distance(*proposed_included, *first_excluded) : true;
    }
    
    template <typename BidirectionalIterator, typename SampleReadMap>
    GenomicRegion
    expand_around_included(BidirectionalIterator first_previous, BidirectionalIterator first_included,
                           BidirectionalIterator first_excluded, BidirectionalIterator last,
                           const SampleReadMap& reads)
    {
        auto last_included = std::prev(first_excluded);
        
        auto leftmost_region  = get_region(*leftmost_overlapped(reads, *first_included, MappableRangeOrder::BidirectionallySorted));
        auto rightmost_region = get_region(*rightmost_overlapped(reads, *last_included, MappableRangeOrder::BidirectionallySorted));
        
        if (count_overlapped(first_previous, first_included, leftmost_region, MappableRangeOrder::BidirectionallySorted) > 0) {
            auto max_left_flank_size = inner_distance(*first_included, *rightmost_mappable(first_previous, first_included));
            leftmost_region = shift(get_region(*first_included), max_left_flank_size);
            
            if (overlaps(*std::prev(first_included), leftmost_region)) {
                // to deal with case where last previous is insertion, otherwise would be included in overlap_range
                leftmost_region = shift(leftmost_region, 1);
            }
        }
        
        if (first_excluded != last && contains(rightmost_region, *first_excluded)) {
            auto max_right_flank_size = inner_distance(*last_included, *first_excluded);
            rightmost_region = shift(get_region(*last_included), max_right_flank_size);
        }
        
        return get_closed(leftmost_region, rightmost_region);
    }
    
} // end namespace detail

enum class IndicatorLimit { SharedWithPreviousRegion, NoLimit };
enum class ExtensionLimit { WithinReadLengthOfFirstIncluded, NoLimit };

/**
 Determines the optimal sub-region given the conditions. Up to 'max_included' non-overlapping variants are included
 in the returned region, of which up to 'max_indicators' may be present in 'the_previous_sub_region'.
 
 The returned region may include less than 'max_included' variants if including more variants would
 cause the region density to increase more than the available number of variants.
 
 The algorithm only uses indicators that may share reads containing variants with the new region (note
 this does not actually look at the read sequence, only the position). Unless 
 'limit_indicators_to_shared_with_previous_region' is false, in which case all indicators are used.
 */
template <typename SampleReadMap, typename Container>
GenomicRegion advance_region(const GenomicRegion& previous_region,
                             const SampleReadMap& reads, const Container& variants,
                             unsigned max_included, unsigned max_indicators,
                             IndicatorLimit indicator_limit=IndicatorLimit::SharedWithPreviousRegion,
                             ExtensionLimit extension_limit=ExtensionLimit::WithinReadLengthOfFirstIncluded)
{
    if (max_included > 0 && max_included <= max_indicators) {
        max_indicators = max_included - 1;
    }
    
    auto candidate_order = MappableRangeOrder::BidirectionallySorted;
    auto read_order      = MappableRangeOrder::BidirectionallySorted;
    
    auto last_variant_it = std::cend(variants);
    
    auto previous_variants = bases(overlap_range(std::cbegin(variants), last_variant_it, previous_region,
                                                 candidate_order));
    auto first_previous_it = previous_variants.begin();
    auto included_it       = previous_variants.end();
    
    if (max_included == 0) {
        return (included_it != last_variant_it) ? get_intervening(previous_region, *included_it) :
                                                    previous_region;
    }
    
    unsigned num_indicators {max_indicators};
    
    if (indicator_limit == IndicatorLimit::SharedWithPreviousRegion) {
        auto first_shared_in_previous_range_it = find_first_shared(reads, first_previous_it,
                                                                   included_it, *included_it,
                                                                   read_order);
        auto num_possible_indicators = static_cast<unsigned>(std::distance(first_shared_in_previous_range_it, included_it));
        
        num_indicators = std::min(num_possible_indicators, num_indicators);
    }
    
    max_included -= num_indicators;
    auto first_included_it = std::prev(included_it, num_indicators);
    
    auto num_remaining_variants = static_cast<unsigned>(std::distance(included_it, last_variant_it));
    unsigned num_excluded_variants {0};
    
    if (extension_limit == ExtensionLimit::WithinReadLengthOfFirstIncluded) {
        auto max_num_variants_within_read_length = static_cast<unsigned>(max_count_if_shared_with_first(reads,
                                                                    included_it, last_variant_it, read_order));
        max_included = std::min({max_included, num_remaining_variants, max_num_variants_within_read_length + 1});
        num_excluded_variants = max_num_variants_within_read_length - max_included;
    } else {
        max_included = std::min(max_included, num_remaining_variants);
    }
    
    auto first_excluded_it = std::next(included_it, max_included);
    
    while (--max_included > 0 &&
           detail::is_optimal_to_extend(first_included_it, std::next(included_it), first_excluded_it,
                                        last_variant_it, reads, max_included + num_excluded_variants)) {
        ++included_it;
    }
    
    std::advance(included_it, count_overlapped(std::next(included_it), last_variant_it,
                                               *rightmost_mappable(first_included_it, std::next(included_it)),
                                               candidate_order));
    
    first_excluded_it = std::next(included_it);
    
    return detail::expand_around_included(first_previous_it, first_included_it, first_excluded_it,
                                          last_variant_it, reads);
}

template <typename SampleReadMap, typename Container>
std::vector<GenomicRegion> cover_region(const GenomicRegion& region,
                                        const SampleReadMap& reads, const Container& variants,
                                        unsigned max_included, unsigned max_indicators,
                                        IndicatorLimit indicator_limit=IndicatorLimit::SharedWithPreviousRegion,
                                        ExtensionLimit extension_limit=ExtensionLimit::WithinReadLengthOfFirstIncluded)
{
    std::vector<GenomicRegion> result {};
    
    auto sub_region = compress_right(region, -size(region));
    
    std::cout << sub_region << std::endl;
    
    while (ends_before(sub_region, region)) {
        sub_region = advance_region(sub_region, reads, variants, max_included, max_indicators,
                                    indicator_limit, extension_limit);
        result.emplace_back(sub_region);
    }
    
    result.shrink_to_fit();
    
    return result;
}
    
} // end namespace Octopus

#endif
