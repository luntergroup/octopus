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
#include <cassert>
#include <limits>

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
        Model model {2, 3}; // TODO: make part of an error model
        
        auto gap_open_penalities = indel_error_model_.calculate_gap_open_penalties(haplotype);
        
        // TODO: we can cache the kmer-hashes to speed this up
        const auto best_mapping_positions = map_query_to_target<5>(read.get_sequence(), haplotype.get_sequence());
        //std::vector<std::size_t> best_mapping_positions {};
        
        auto max_log_probability = std::numeric_limits<double>::lowest();
        
        for (const auto position : best_mapping_positions) {
            auto cur = compute_log_conditional_probability(haplotype.get_sequence(), read.get_sequence(),
                                                           read.get_qualities(), gap_open_penalities,
                                                           position, model);
            if (cur > max_log_probability) max_log_probability = cur;
        }
        
        assert(contains(haplotype, read));
        
        const auto original_mapping_position = begin_distance(read, haplotype);
        
        assert(std::is_sorted(std::cbegin(best_mapping_positions), std::cend(best_mapping_positions)));
        
        if (!std::binary_search(std::cbegin(best_mapping_positions),
                                std::cend(best_mapping_positions),
                                original_mapping_position) && max_log_probability < 0) {
            auto cur = compute_log_conditional_probability(haplotype.get_sequence(), read.get_sequence(),
                                                           read.get_qualities(), gap_open_penalities,
                                                           original_mapping_position, model);
            
            if (cur > max_log_probability) max_log_probability = cur;
        }
        
        return max_log_probability;
    }
} // namespace Octopus
