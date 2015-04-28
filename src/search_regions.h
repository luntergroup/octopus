//
//  search_regions.h
//  Octopus
//
//  Created by Daniel Cooke on 27/04/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_search_regions_h
#define Octopus_search_regions_h

#include <iterator>  // std::distance, std::next, std::prev, std::cbegin, std::cend
#include <algorithm> // std::min
#include <cmath>     // std::abs

#include "genomic_region.h"
#include "region_utils.h"
#include "map_utils.h"

#include <iostream> // TEST
using std::cout;    // TEST
using std::endl;    // TEST

/**
 A region is 'dense' if the number of variants within a read length of the first variant in
 'the_variants' that overlaps 'a_region', is greater than max_variants
 */
template <typename SampleReadMap, typename Variants>
bool is_dense_region(const GenomicRegion& a_region, const SampleReadMap& the_reads,
                     const Variants& the_variants, unsigned max_variants)
{
    auto overlapped = overlap_range(std::cbegin(the_variants), std::cend(the_variants), a_region);
    
    auto max_num_variants_within_read_length = max_count_if_shared_with_first(the_reads, overlapped.first,
                                                                              std::cend(the_variants));
    
    return max_num_variants_within_read_length > max_variants;
}

namespace detail
{
    template <typename BidirectionalIterator, typename SampleReadMap>
    bool is_optimal_to_include(BidirectionalIterator first_included, BidirectionalIterator current,
                               BidirectionalIterator first_excluded, BidirectionalIterator last,
                               const SampleReadMap& the_reads, unsigned max_density_increase)
    {
        if (current == first_included) return true;
        
        bool increases_density {max_count_if_shared_with_first(the_reads, std::next(current), last)
                                        > max_density_increase};
        
        return (increases_density) ? outer_distance(*std::prev(current), *current)
                                        <= outer_distance(*current, *first_excluded) : true;
    }
    
    template <typename BidirectionalIterator, typename SampleReadMap>
    GenomicRegion
    get_optimal_region_around_variants(BidirectionalIterator first_previous, BidirectionalIterator first_included,
                                       BidirectionalIterator first_excluded, BidirectionalIterator last,
                                       const SampleReadMap& the_reads)
    {
        auto last_included = std::prev(first_excluded);
        
        auto leftmost_region  = leftmost_overlapped(the_reads, *first_included)->get_region();
        auto rightmost_region = rightmost_overlapped(the_reads, *last_included)->get_region();
        
        if (count_overlapped(first_previous, first_included, leftmost_region) > 0) {
            auto max_left_flank_size = inner_distance(*rightmost_mappable(first_previous, first_included),
                                                      *first_included);
            leftmost_region = shift(first_included->get_region(), -max_left_flank_size);
        }
        
        if (first_excluded != last && contains(rightmost_region, *first_excluded)) {
            auto max_right_flank_size = inner_distance(*last_included, *first_excluded);
            rightmost_region = shift(last_included->get_region(), max_right_flank_size);
        }
        
        return get_closed(leftmost_region, rightmost_region);
    }
    
} // end namespace detail

template <typename SampleReadMap, typename Variants>
GenomicRegion next_sub_region(const GenomicRegion& the_search_region, const GenomicRegion& the_previous_sub_region,
                              const SampleReadMap& the_reads, const Variants& the_variants,
                              unsigned max_variants, unsigned max_indicators)
{
    if (max_variants > 0 && max_variants <= max_indicators) {
        max_indicators = max_variants - 1;
    }
    
    auto last_variant_it  = std::cend(the_variants);
    
    auto previous_variant_sub_range = overlap_range(std::cbegin(the_variants), last_variant_it, the_previous_sub_region);
    auto included_it = previous_variant_sub_range.second;
    
    auto first_shared_in_previous_range_it = find_first_shared(the_reads, previous_variant_sub_range.first,
                                                               previous_variant_sub_range.second, *included_it);
    auto num_possible_indicators = static_cast<unsigned>(std::distance(first_shared_in_previous_range_it, included_it));
    unsigned num_indicators = std::min(num_possible_indicators, max_indicators);
    
    max_variants -= num_indicators;
    auto first_included_it = std::prev(included_it, num_indicators);
    
    auto num_remaining_variants = static_cast<unsigned>(std::distance(included_it, last_variant_it));
    if (num_remaining_variants < max_variants) {
        max_variants -= num_remaining_variants;
    }
    
    auto max_num_variants_within_read_length = static_cast<unsigned>(max_count_if_shared_with_first(the_reads,
                                                                    included_it, last_variant_it));
    max_variants = std::min(max_variants, max_num_variants_within_read_length + 1);
    
    unsigned num_excluded_variants = max_num_variants_within_read_length - max_variants;
    auto first_excluded_it = std::next(included_it, max_variants);
    
    // Any 'excluded' variants overlapping the last included variant will be implicitly included in
    // the final region, so they are added explictly here
    while (first_excluded_it != last_variant_it && overlaps(*std::prev(first_excluded_it), *first_excluded_it)) {
        ++first_excluded_it;
        ++max_variants;
    }
    
    while (included_it != first_excluded_it &&
           detail::is_optimal_to_include(first_included_it, included_it, first_excluded_it,
                                         last_variant_it, the_reads,
                                         max_variants + num_excluded_variants)) {
        ++included_it;
    }
    
    return detail::get_optimal_region_around_variants(previous_variant_sub_range.first, first_included_it,
                                                      included_it, last_variant_it, the_reads);
}

#endif
