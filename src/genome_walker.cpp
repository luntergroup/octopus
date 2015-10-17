//
//  genome_walker.cpp
//  Octopus
//
//  Created by Daniel Cooke on 17/10/2015.
//  Copyright © 2015 Oxford University. All rights reserved.
//

#include "genome_walker.hpp"

#include <iterator>  // std::distance, std::next, std::prev, std::advance, std::cbegin, std::cend
#include <algorithm> // std::min
#include <cmath>     // std::abs

#include "genomic_region.hpp"
#include "mappable.hpp"
#include "aligned_read.hpp"
#include "variant.hpp"
#include "mappable_algorithms.hpp"
#include "mappable_map.hpp"
#include "map_utils.hpp"

#include <iostream> // DEBUG

GenomeWalker::GenomeWalker(unsigned max_indicators, unsigned max_included,
                           IndicatorLimit indicator_limit, ExtensionLimit extension_limit)
:
max_indicators_ {(max_included > 0 && max_included <= max_indicators) ? max_included - 1 : max_indicators},
max_included_ {max_included},
indicator_limit_ {indicator_limit},
extension_limit_ {extension_limit}
{}

GenomicRegion GenomeWalker::start_walk(const ContigNameType& contig, const ReadMap& reads,
                                       const Candidates& candidates)
{
    return walk(GenomicRegion {contig, 0, 0}, reads, candidates);
}

GenomicRegion GenomeWalker::continue_walk(const GenomicRegion& previous_region, const ReadMap& reads,
                                          const Candidates& candidates)
{
    return walk(previous_region, reads, candidates);
}

template <typename BidirectionalIterator, typename SampleReadMap>
bool is_optimal_to_extend(BidirectionalIterator first_included, BidirectionalIterator proposed_included,
                          BidirectionalIterator first_excluded, BidirectionalIterator last,
                          const SampleReadMap& reads, unsigned max_density_increase)
{
    if (proposed_included == last) return false;
    if (first_excluded == last) return true;
    
    bool increases_density {max_count_if_shared_with_first(reads, proposed_included, last)
        >= max_density_increase};
    
    return !increases_density || inner_distance(*std::prev(proposed_included), *proposed_included) <= inner_distance(*proposed_included, *first_excluded);
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

GenomicRegion GenomeWalker::walk(const GenomicRegion& previous_region, const ReadMap& reads,
                                 const Candidates& candidates)
{
    using std::cbegin; using std::cend; using std::next; using std::prev; using std::min;
    using std::distance; using std::advance;
    
    //std::cout << "walking from " << previous_region << std::endl;
    
    auto last_variant_it = cend(candidates);
    
    auto previous_candidates = bases(candidates.overlap_range(previous_region));
    
    auto first_previous_it = cbegin(previous_candidates);
    auto included_it       = cend(previous_candidates);
    
    if (included_it == last_variant_it) return get_tail(previous_region, 0);
    
    if (max_included_ == 0) {
        if (included_it != last_variant_it) {
            return get_intervening(previous_region, *included_it);
        } else {
            return previous_region;
        }
    }
    
    unsigned num_indicators {max_indicators_};
    
    if (indicator_limit_ == IndicatorLimit::SharedWithPreviousRegion) {
        auto first_shared_in_previous_range_it = find_first_shared(reads, first_previous_it,
                                                                   included_it, *included_it);
        auto num_possible_indicators = static_cast<unsigned>(distance(first_shared_in_previous_range_it, included_it));
        
        num_indicators = min(num_possible_indicators, num_indicators);
    }
    
    auto num_included = max_included_ - num_indicators;
    
    auto first_included_it = prev(included_it, num_indicators);
    
    auto num_remaining_candidates = static_cast<unsigned>(distance(included_it, last_variant_it));
    unsigned num_excluded_candidates {0};
    
    if (extension_limit_ == ExtensionLimit::WithinReadLengthOfFirstIncluded) {
        auto max_num_candidates_within_read_length = static_cast<unsigned>(max_count_if_shared_with_first(reads,
                                                                                                        included_it, last_variant_it));
        num_included = min({num_included, num_remaining_candidates, max_num_candidates_within_read_length + 1});
        num_excluded_candidates = max_num_candidates_within_read_length - num_included;
    } else {
        num_included = min(num_included, num_remaining_candidates);
    }
    
    auto first_excluded_it = next(included_it, num_included);
    
    while (--num_included > 0 &&
           is_optimal_to_extend(first_included_it, next(included_it), first_excluded_it,
                                last_variant_it, reads, num_included + num_excluded_candidates)) {
               ++included_it;
           }
    
    advance(included_it, candidates.count_overlapped(*rightmost_mappable(first_included_it, next(included_it))) - 1);
    
    first_excluded_it = next(included_it);
    
    auto lhs_read = *leftmost_overlapped(reads, *first_included_it);
    auto rhs_read = *rightmost_overlapped(reads, *included_it);
    
    return get_encompassing(lhs_read, rhs_read);
}
