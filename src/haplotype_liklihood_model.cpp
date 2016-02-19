//
//  haplotype_liklihood_model.cpp
//  Octopus
//
//  Created by Daniel Cooke on 25/09/2015.
//  Copyright © 2015 Oxford University. All rights reserved.
//

#include "haplotype_liklihood_model.hpp"

#include <cmath>
#include <algorithm>

#include "mappable.hpp"
#include "aligned_read.hpp"
#include "kmer_mapping.hpp"
#include "pair_hmm.hpp"

#include <iostream> // TEST
#include <chrono>   // TEST

namespace Octopus
{
    // public methods
    
    double HaplotypeLikelihoodModel::log_probability(const AlignedRead& read, const Haplotype& haplotype)
    {
        Model model {2, 3};
        
        auto gap_open_penalities = indel_error_model_.calculate_gap_open_penalties(haplotype);
        
        //const auto offset_hints = map_query_to_target<5>(read.get_sequence(), haplotype.get_sequence());
        
        auto offset_hint = begins_before(haplotype, read) ? begin_distance(read, haplotype) : 0;
        
//        if (!offset_hints.empty() && offset_hints.front() != offset_hint) {
//            offset_hint = offset_hints.front();
//        }
        
        return compute_log_conditional_probability(haplotype.get_sequence(), read.get_sequence(),
                                                   read.get_qualities(), gap_open_penalities,
                                                   offset_hint, model);
    }
} // namespace Octopus
