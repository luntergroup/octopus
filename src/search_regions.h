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
using std::cout;    // TEST
using std::endl;    // TEST

namespace Octopus
{

/**
 A region is 'dense' if the number of variants within a read length of the first variant in
 'the_variants' that overlaps 'a_region', is greater than max_variants
 */
template <typename SampleReadMap, typename Container>
bool is_dense_region(const GenomicRegion& a_region, const SampleReadMap& the_reads,
                     const Container& the_variants, unsigned max_variants)
{
    auto overlapped = overlap_range(std::cbegin(the_variants), std::cend(the_variants), a_region);
    
    auto max_num_variants_within_read_length = max_count_if_shared_with_first(the_reads, overlapped.first,
                                                                              std::cend(the_variants));
    
    return max_num_variants_within_read_length > max_variants;
}

namespace detail
{
    template <typename BidirectionalIterator, typename SampleReadMap>
    bool is_optimal_to_extend(BidirectionalIterator first_included, BidirectionalIterator proposed_included,
                              BidirectionalIterator first_excluded, BidirectionalIterator last,
                              const SampleReadMap& the_reads, unsigned max_density_increase)
    {
        if (proposed_included == last) return false;
        if (first_excluded == last) return true;
        
        bool increases_density {max_count_if_shared_with_first(the_reads, proposed_included, last)
                                        >= max_density_increase};
        
        return (increases_density) ? inner_distance(*std::prev(proposed_included), *proposed_included)
                                        <= inner_distance(*proposed_included, *first_excluded) : true;
    }
    
    template <typename BidirectionalIterator, typename SampleReadMap>
    GenomicRegion
    get_optimal_region_around_included(BidirectionalIterator first_previous, BidirectionalIterator first_included,
                                       BidirectionalIterator first_excluded, BidirectionalIterator last,
                                       const SampleReadMap& the_reads)
    {
        auto last_included = std::prev(first_excluded);
        
        auto leftmost_region  = leftmost_overlapped(the_reads, *first_included)->get_region();
        auto rightmost_region = rightmost_overlapped(the_reads, *last_included)->get_region();

        if (count_overlapped(first_previous, first_included, leftmost_region) > 0) {
            auto max_left_flank_size = inner_distance(*first_included, *rightmost_mappable(first_previous, first_included));
            leftmost_region = shift(first_included->get_region(), max_left_flank_size);
            
            if (overlaps(*std::prev(first_included), leftmost_region)) {
                // to deal with case where last previous is insertion, otherwise would be included in overlap_range
                leftmost_region = shift(leftmost_region, 1);
            }
        }
        
        if (first_excluded != last && contains(rightmost_region, *first_excluded)) {
            auto max_right_flank_size = inner_distance(*last_included, *first_excluded);
            rightmost_region = shift(last_included->get_region(), max_right_flank_size);
        }
        
        return get_closed(leftmost_region, rightmost_region);
    }
    
} // end namespace detail

template <typename SampleReadMap, typename Container>
GenomicRegion next_sub_region(const GenomicRegion& the_previous_sub_region,
                              const SampleReadMap& the_reads, const Container& the_variants,
                              unsigned max_included, unsigned max_indicators)
{
    if (max_included > 0 && max_included <= max_indicators) {
        max_indicators = max_included - 1;
    }
    
    auto last_variant_it  = std::cend(the_variants);
    
    auto previous_variant_sub_range = overlap_range(std::cbegin(the_variants), last_variant_it, the_previous_sub_region);
    auto included_it = previous_variant_sub_range.second;
    
    if (max_included == 0) {
        return (included_it != last_variant_it) ? get_intervening(the_previous_sub_region, *included_it) :
                                                    the_previous_sub_region;
    }
    
    auto first_shared_in_previous_range_it = find_first_shared(the_reads, previous_variant_sub_range.first,
                                                               previous_variant_sub_range.second, *included_it);
    auto num_possible_indicators = static_cast<unsigned>(std::distance(first_shared_in_previous_range_it, included_it));
    
    unsigned num_indicators = std::min(num_possible_indicators, max_indicators);
    
    max_included -= num_indicators;
    auto first_included_it = std::prev(included_it, num_indicators);
    
    auto num_remaining_variants = static_cast<unsigned>(std::distance(included_it, last_variant_it));
    auto max_num_variants_within_read_length = static_cast<unsigned>(max_count_if_shared_with_first(the_reads,
                                                                    included_it, last_variant_it));
    
    max_included = std::min({max_included, num_remaining_variants, max_num_variants_within_read_length + 1});
    
    unsigned num_excluded_variants = max_num_variants_within_read_length - max_included;
    auto first_excluded_it = std::next(included_it, max_included);
    
    while (--max_included > 0 &&
           detail::is_optimal_to_extend(first_included_it, std::next(included_it), first_excluded_it,
                                        last_variant_it, the_reads, max_included + num_excluded_variants)) {
        ++included_it;
    }
    
    std::advance(included_it, count_overlapped(std::next(included_it), last_variant_it,
                                               *rightmost_mappable(first_included_it, std::next(included_it))));
    
    first_excluded_it = std::next(included_it);
    
    return detail::get_optimal_region_around_included(previous_variant_sub_range.first, first_included_it,
                                                      first_excluded_it, last_variant_it, the_reads);
}
    
} // end namespace Octopus

#endif
