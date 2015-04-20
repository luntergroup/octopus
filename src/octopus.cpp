//
//  octopus.cpp
//  Octopus
//
//  Created by Daniel Cooke on 08/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "octopus.h"

#include <iterator>
#include <algorithm>

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

GenomicRegion next_sub_region(const GenomicRegion& the_search_region, const GenomicRegion& previous_sub_region,
                              const ReadManager::SampleReadMap& the_reads, const Variants& the_candidates,
                              unsigned max_variants_in_region, unsigned max_region_size, unsigned max_region_overlap)
{
    auto previous_sub_range = overlap_range(std::cbegin(the_candidates), std::cend(the_candidates), previous_sub_region);
    
    auto included_candidate_it = previous_sub_range.second;
    
    if (!empty(previous_sub_region)) {
        auto num_samples_with_shared_read_candidates = std::count_if(std::cbegin(the_reads), std::cend(the_reads),
                               [previous_sub_range, included_candidate_it] (const auto& sample_reads) {
                                   return has_shared(std::cbegin(sample_reads.second),
                                                     std::cend(sample_reads.second),
                                                     *previous_sub_range.second, *included_candidate_it);
                               });
        
        if (num_samples_with_shared_read_candidates > 0) {
            //--included_candidate_it;
        }
    }
    
    auto first_included_it = included_candidate_it;
    bool shares_reads_with_next_candidate {true};
    auto next_candidate_it = included_candidate_it;
    
    --max_variants_in_region;
    
    while (max_variants_in_region > 0 && shares_reads_with_next_candidate) {
        
        cout << *included_candidate_it << endl;
        
        next_candidate_it = std::next(included_candidate_it);
        
        if (next_candidate_it == std::cend(the_candidates)) {
            break;
        }
        
        shares_reads_with_next_candidate = std::any_of(std::cbegin(the_reads), std::cend(the_reads),
                      [included_candidate_it, next_candidate_it] (const auto& sample_reads) {
                          return has_shared(std::cbegin(sample_reads.second), std::cend(sample_reads.second),
                                            *included_candidate_it, *next_candidate_it);
                      });
        
        included_candidate_it = next_candidate_it;
        --max_variants_in_region;
    }
    
    std::vector<GenomicRegion> right_flank_read_regions {};
    
    std::transform(std::cbegin(the_reads), std::cend(the_reads), std::back_inserter(right_flank_read_regions),
                   [included_candidate_it] (const auto& sample_reads) {
                       return std::prev(overlap_range(std::cbegin(sample_reads.second), std::cend(sample_reads.second),
                                                      *included_candidate_it).second)->get_region();
                   });
    
    std::sort(right_flank_read_regions.begin(), right_flank_read_regions.end());
    
    
    
    return get_encompassing_region(*first_included_it, *included_candidate_it);
}
