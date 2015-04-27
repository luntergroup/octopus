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

/**
 This algorithm finds the 'optimal' next sub-region with the search region 'the_search_region',
 where the optimal sub-region is the one that will minimise the number of sub-regions within
 the search region given the contraints 'max_variants' and 'max_indicators'.
 
 The algorithm may return a region containing less than 'max_variants' if including more variants
 into the region would cause a non-dense region to become dense.
 
 The condition 'max_indicators' bounds how many variants in the previous region can be included
 in the new region if they share overlapping reads. This is crucial for haplotype calling in dense 
 regions where the haplotype is split into multiple sub-regions.
 */
template <typename SampleReadMap, typename Variants>
GenomicRegion next_sub_region(const GenomicRegion& the_search_region,
                              const GenomicRegion& the_previous_sub_region,
                              const SampleReadMap& the_reads, const Variants& the_variants,
                              unsigned max_variants, unsigned max_indicators)
{
    if (max_variants <= max_indicators) {
        max_indicators = max_variants - 1;
    }
    
    auto previous_variant_sub_range = overlap_range(std::cbegin(the_variants), std::cend(the_variants),
                                                    the_previous_sub_region);
    
    auto included_it = previous_variant_sub_range.second;
    
    auto first_shared_in_previous_range_it = find_first_shared(the_reads, previous_variant_sub_range.first,
                                                               previous_variant_sub_range.second, *included_it);
    
    auto num_possible_indicators = std::distance(first_shared_in_previous_range_it, included_it);
    auto num_indicators = std::min(static_cast<unsigned>(num_possible_indicators), max_indicators);
    
    max_variants -= num_indicators;
    
    auto max_num_variants_within_read_length = max_count_if_shared_with_first(the_reads, included_it,
                                                                              std::cend(the_variants));
    
    max_variants = std::min(max_variants, static_cast<unsigned>(max_num_variants_within_read_length));
    
    auto num_excluded_variants = max_num_variants_within_read_length - max_variants;
    
    auto first_new_included_it = included_it;
    
    while (max_variants > 1 && max_count_if_shared_with_first(the_reads, std::next(included_it),
                                    std::cend(the_variants)) < max_variants + num_excluded_variants) {
        --max_variants;
        ++included_it;
    }
    
    cout << *first_new_included_it << endl;
    cout << *included_it << endl;
    
    auto leftmost_region  = leftmost_overlapping(the_reads, *first_new_included_it)->get_region();
    auto rightmost_region = rightmost_overlapping(the_reads, *included_it)->get_region();
    
    if (num_indicators == 0 && num_possible_indicators > 0) {
        auto trim_size = -distance(leftmost_region, *std::prev(previous_variant_sub_range.second));
        
        if (std::abs(trim_size) <= size(leftmost_region)) {
            leftmost_region = compress_left(leftmost_region, trim_size);
        } else {
            auto max_trim = std::min(size(leftmost_region), static_cast<GenomicRegion::SizeType>(trim_size));
            rightmost_region = shift(compress_left(leftmost_region, max_trim), trim_size - max_trim);
        }
    }
    
    if (max_variants > 0) {
        auto trim_size = distance(*std::next(included_it), rightmost_region);
        
        if (std::abs(trim_size) <= size(rightmost_region)) {
            rightmost_region = compress_right(rightmost_region, trim_size);
        } else {
            auto max_trim = std::min(size(rightmost_region), static_cast<GenomicRegion::SizeType>(trim_size));
            rightmost_region = shift(compress_right(rightmost_region, max_trim), trim_size - max_trim);
        }
    }
    
    return get_encompassing_region(leftmost_region, rightmost_region);
}

#endif
