//
//  octopus.cpp
//  Octopus
//
//  Created by Daniel Cooke on 08/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "octopus.h"

#include <iterator>  // std::distance, std::cbegin etc
#include <algorithm> // std::min

#include "genomic_region.h"
#include "aligned_read.h"
#include "variant.h"
#include "region_utils.h"

#include <iostream> // TEST
using std::cout; // TEST
using std::endl; // TEST

void run_octopus()
{
    
}

GenomicRegion next_sub_region(const GenomicRegion& the_search_region, const GenomicRegion& the_previous_sub_region,
                              const ReadManager::SampleReadMap& the_reads, const Variants& the_candidates,
                              unsigned max_variants, unsigned max_indicators)
{
    if (max_variants <= max_indicators) {
        max_indicators = max_indicators - 1;
    }
    
    auto previous_sub_range = overlap_range(std::cbegin(the_candidates), std::cend(the_candidates),
                                            the_previous_sub_region);
    
    auto included_it = previous_sub_range.second;
    
    auto first_included_it = included_it;
    
    auto first_shared_in_previous_range_it = find_first_shared(the_reads, previous_sub_range.first,
                                                               previous_sub_range.second, *included_it);
    
    auto num_possible_indicators = std::distance(first_shared_in_previous_range_it, included_it);
    
    while (num_possible_indicators > 0 && max_indicators > 0) {
        --included_it;
        --num_possible_indicators;
        --max_indicators;
        --max_variants;
    }
    
    auto last_candidate_it = std::cend(the_candidates);
    
    auto max_num_variants_within_read_length = max_count_if_shared_with_first(the_reads, included_it,
                                                                              last_candidate_it);
    
    max_variants = std::min(max_variants, static_cast<unsigned>(max_num_variants_within_read_length));
    
    auto num_excluded_variants = max_num_variants_within_read_length - max_variants;
    
    while (max_variants > 0 && max_count_if_shared_with_first(the_reads, std::next(included_it),
                                        last_candidate_it) < max_variants + num_excluded_variants) {
        --max_variants;
        ++included_it;
    }
    
    auto leftmost_region  = leftmost_overlapping(the_reads, *first_included_it);
    auto rightmost_region = rightmost_overlapping(the_reads, *included_it);
    
    if (!is_after(leftmost_region, the_previous_sub_region)) {
        leftmost_region = shift(leftmost_region, -distance(the_previous_sub_region, leftmost_region));
    }
    
    if (max_variants > 0) {
        rightmost_region = shift(rightmost_region, distance(rightmost_region, *std::next(included_it)));
    }
    
    return get_encompassing_region(leftmost_region, rightmost_region);
}
