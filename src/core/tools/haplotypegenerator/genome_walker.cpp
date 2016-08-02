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

#include <basics/genomic_region.hpp>
#include <interfaces/mappable.hpp>
#include <basics/aligned_read.hpp>
#include <core/types/variant.hpp>
#include <utils/mappable_algorithms.hpp>
#include <containers/mappable_map.hpp>

#include <iostream> // DEBUG

namespace octopus { namespace coretools {

GenomeWalker::GenomeWalker(unsigned max_included,
                           IndicatorPolicy indicator_policy,
                           ExtensionPolicy extension_policy)
:
max_included_ {max_included},
indicator_policy_ {indicator_policy},
extension_policy_ {extension_policy}
{}

GenomicRegion GenomeWalker::walk(const ContigName& contig, const ReadMap& reads,
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
    
    if (included_itr == last_candidate_itr) {
        return shift(tail_region(previous_region), 2);
    }
    
    if (max_included_ == 0) {
        if (included_itr != last_candidate_itr) {
            return intervening_region(previous_region, *included_itr);
        } else {
            return previous_region;
        }
    }
    
    unsigned num_indicators {0};
    
    switch (indicator_policy_) {
        case IndicatorPolicy::IncludeNone: break;
        case IndicatorPolicy::IncludeIfSharedWithNovelRegion:
        {
            if (distance(first_previous_itr, included_itr) > 0) {
                const auto it = find_first_shared(reads, first_previous_itr, included_itr, *included_itr);
                
                num_indicators = static_cast<unsigned>(distance(it, included_itr));
            }
            break;
        }
        case IndicatorPolicy::IncludeIfLinkableToNovelRegion:
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
        case IndicatorPolicy::IncludeAll:
            num_indicators = static_cast<unsigned>(distance(first_previous_itr, included_itr));
            break;
    }
    
    auto first_included_itr = prev(included_itr, num_indicators);
    
    auto num_remaining_candidates = static_cast<unsigned>(distance(included_itr, last_candidate_itr));
    
    unsigned num_excluded_candidates {0};
    
    auto num_included = max_included_;
    
    if (extension_policy_ == ExtensionPolicy::IncludeIfWithinReadLengthOfFirstIncluded) {
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
               if (extension_policy_ == ExtensionPolicy::IncludeIfSharedWithFrontier
                   && !has_shared(reads, *included_itr, *next(included_itr))) {
                   break;
               }
               ++included_itr;
           }
    
    // first_excluded_itr = find_first_after(included_itr, last_candidate_itr, rightmost
    
    first_excluded_itr = next(included_itr);
    
    const auto& rightmost = *rightmost_mappable(first_included_itr, first_excluded_itr);
    
    auto num_remaining = static_cast<size_t>(distance(first_excluded_itr, last_candidate_itr));
    
    auto num_overlapped = candidates.count_overlapped(first_excluded_itr, last_candidate_itr, rightmost);
    
    advance(included_itr, min(num_remaining, num_overlapped));
    
    return encompassing_region(first_included_itr, included_itr);
}

} // namespace coretools
} // namespace octopus
