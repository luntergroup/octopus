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

namespace details {
    
    template <typename BidirectionalIterator, typename SampleReadMap>
    bool is_optimal_to_include(BidirectionalIterator current, BidirectionalIterator first_out_of_bounds,
                               BidirectionalIterator last, const SampleReadMap& the_reads,
                               unsigned max_density_increase)
    {
        bool increases_density {max_count_if_shared_with_first(the_reads, std::next(current), last)
                                        > max_density_increase};
        
        cout << *current << endl;
        cout << distance(*std::prev(current), *current) << endl;
        cout << distance(*current, *first_out_of_bounds) << endl;
        
        return (increases_density) ? distance(*std::prev(current), *current)
                                        <= distance(*current, *first_out_of_bounds) : true;
    }
    
    template <typename BidirectionalIterator, typename SampleReadMap, typename Variants>
    GenomicRegion
    get_optimal_region_around_variants(BidirectionalIterator first, BidirectionalIterator last,
                                       const SampleReadMap& the_reads, const Variants& the_variants,
                                       unsigned num_used_indicators, unsigned num_unused_variants)
    {
        auto leftmost_region  = leftmost_overlapping(the_reads, *first)->get_region();
        auto rightmost_region = rightmost_overlapping(the_reads, *last)->get_region();
        
        //cout << rightmost_region << endl;
        
        if (num_used_indicators == 0) {
            auto trim_size = distance(leftmost_region, *std::prev(first));
            
            if (std::abs(trim_size) <= size(leftmost_region)) {
                leftmost_region = compress_left(leftmost_region, trim_size);
            } else {
                auto max_trim = std::min(size(leftmost_region), static_cast<GenomicRegion::SizeType>(trim_size));
                rightmost_region = shift(compress_left(leftmost_region, max_trim), trim_size - max_trim);
            }
        }
        
        if (num_unused_variants > 0) {
            auto trim_size = -distance(*std::next(last), rightmost_region);
            
            if (std::abs(trim_size) <= size(rightmost_region)) {
                rightmost_region = compress_right(rightmost_region, trim_size);
            } else {
                auto max_trim = std::min(size(rightmost_region), static_cast<GenomicRegion::SizeType>(trim_size));
                rightmost_region = shift(compress_right(rightmost_region, max_trim), trim_size - max_trim);
            }
        }
        
        return get_encompassing_region(leftmost_region, rightmost_region);
    }
    
} // end namespace details

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
    
    auto last_variant_it = std::cend(the_variants);
    auto previous_variant_sub_range = overlap_range(std::cbegin(the_variants), last_variant_it,
                                                    the_previous_sub_region);
    
    auto included_it = previous_variant_sub_range.second;
    
    auto first_shared_in_previous_range_it = find_first_shared(the_reads, previous_variant_sub_range.first,
                                                               previous_variant_sub_range.second, *included_it);
    auto num_possible_indicators = static_cast<unsigned>(std::distance(first_shared_in_previous_range_it, included_it));
    unsigned num_indicators = std::min(num_possible_indicators, max_indicators);
    
    max_variants -= num_indicators;
    
    auto num_remaining_variants = static_cast<unsigned>(std::distance(included_it, last_variant_it));
    if (num_remaining_variants < max_variants) {
        max_variants -= num_remaining_variants;
    }
    
    auto max_num_variants_within_read_length = static_cast<unsigned>(max_count_if_shared_with_first(the_reads,
                                                                    included_it, last_variant_it));
    max_variants = std::min(max_variants, max_num_variants_within_read_length);
    
    unsigned num_excluded_variants = max_num_variants_within_read_length - max_variants;
    auto first_out_of_bounds_it    = std::next(included_it, max_variants + 1);
    auto first_new_included_it     = included_it;
    
    while (--max_variants > 0 && details::is_optimal_to_include(included_it, first_out_of_bounds_it,
                                last_variant_it, the_reads, max_variants + num_excluded_variants)) {
        ++included_it;
    }
    
//    cout << the_previous_sub_region << endl;
//    cout << *first_new_included_it << endl;
//    cout << *included_it << endl;
    
    return details::get_optimal_region_around_variants(first_new_included_it, included_it, the_reads,
                                                       the_variants, num_indicators, max_variants);
}

#endif
