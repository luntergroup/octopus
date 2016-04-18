//
//  genome_walker.cpp
//  Octopus
//
//  Created by Daniel Cooke on 17/10/2015.
//  Copyright © 2015 Oxford University. All rights reserved.
//

#include "genome_walker.hpp"

#include <iterator>
#include <algorithm>
#include <cmath>

#include "genomic_region.hpp"
#include "mappable.hpp"
#include "aligned_read.hpp"
#include "variant.hpp"
#include "mappable_algorithms.hpp"
#include "mappable_map.hpp"

#include <iostream> // DEBUG

namespace Octopus
{
GenomeWalker::GenomeWalker(unsigned max_included, IndicatorLimit indicator_limit,
                           ExtensionLimit extension_limit)
:
max_included_ {max_included},
indicator_limit_ {indicator_limit},
extension_limit_ {extension_limit}
{}

GenomicRegion GenomeWalker::walk(const ContigNameType& contig, const ReadMap& reads,
                                 const Candidates& candidates) const
{
    return walk(GenomicRegion {contig, 0, 0}, reads, candidates);
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

GenomicRegion GenomeWalker::walk(const GenomicRegion& previous_region, const ReadMap& reads,
                                 const Candidates& candidates) const
{
    using std::cbegin; using std::cend; using std::next; using std::prev; using std::min;
    using std::distance; using std::advance;
    
    auto last_candidate_itr = cend(candidates);
    
    auto previous_candidates = bases(overlap_range(candidates, previous_region));
    
    auto first_previous_itr = cbegin(previous_candidates);
    auto included_itr       = cend(previous_candidates);
    
    if (included_itr == last_candidate_itr) return tail_region(previous_region);
    
    if (max_included_ == 0) {
        if (included_itr != last_candidate_itr) {
            return intervening_region(previous_region, *included_itr);
        } else {
            return previous_region;
        }
    }
    
    unsigned num_indicators {0};
    
    if (indicator_limit_ != IndicatorLimit::None) {
        num_indicators = static_cast<unsigned>(distance(first_previous_itr, included_itr));
    }
    
    if (num_indicators > 0 && indicator_limit_ == IndicatorLimit::SharedWithPreviousRegion) {
        auto it = find_first_shared(reads, first_previous_itr, included_itr, *included_itr);
        auto max_possible_indicators = static_cast<unsigned>(distance(it, included_itr));
        num_indicators = min(max_possible_indicators, num_indicators);
    }
    
    auto num_included = max_included_;
    
    auto first_included_itr = prev(included_itr, num_indicators);
    
    auto num_remaining_candidates = static_cast<unsigned>(distance(included_itr, last_candidate_itr));
    
    unsigned num_excluded_candidates {0};
    
    if (extension_limit_ == ExtensionLimit::WithinReadLengthOfFirstIncluded) {
        auto max_candidates_within_read_length = static_cast<unsigned>(max_count_if_shared_with_first(reads, first_included_itr, last_candidate_itr));
        num_included = min({num_included, num_remaining_candidates, max_candidates_within_read_length + 1});
        num_excluded_candidates = max_candidates_within_read_length - num_included;
    } else {
        num_included = min(num_included, num_remaining_candidates);
    }
    
    auto first_excluded_itr = next(included_itr, num_included);
    
    while (--num_included > 0 &&
           is_optimal_to_extend(first_included_itr, next(included_itr), first_excluded_itr,
                                last_candidate_itr, reads, num_included + num_excluded_candidates)) {
               if (extension_limit_ == ExtensionLimit::SharedWithFrontier
                   && !has_shared(reads, *included_itr, *next(included_itr))) {
                   break;
               }
               ++included_itr;
           }
    
    first_excluded_itr = next(included_itr);
    
    const auto& rightmost = *rightmost_mappable(first_included_itr, first_excluded_itr);
    
    auto num_remaining = static_cast<size_t>(distance(first_excluded_itr, last_candidate_itr));
    
    auto num_overlapped = candidates.count_overlapped(first_excluded_itr, last_candidate_itr, rightmost);
    
    advance(included_itr, min(num_remaining, num_overlapped));
    
    first_excluded_itr = next(included_itr);
    
    return encompassing_region(first_included_itr, included_itr);
}

} // namespace Octopus
