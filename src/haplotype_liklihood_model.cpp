//
//  haplotype_liklihood_model.cpp
//  Octopus
//
//  Created by Daniel Cooke on 25/09/2015.
//  Copyright © 2015 Oxford University. All rights reserved.
//

#include "haplotype_liklihood_model.hpp"

#include <utility>
#include <cmath>
#include <algorithm>
#include <cassert>
#include <limits>

#include "mappable.hpp"
#include "aligned_read.hpp"
#include "pair_hmm.hpp"

#include <iostream> // TEST
#include <chrono>   // TEST

namespace Octopus
{
    // public methods
    
    HaplotypeLikelihoodModel::HaplotypeLikelihoodModel(KmerMapper mapper)
    :
    mapper_ {std::move(mapper)}
    {}
    
    namespace
    {
        bool contains(const std::vector<std::size_t>& mapped_positions,
                      const std::size_t original_mapping_position)
        {
            return std::binary_search(std::cbegin(mapped_positions), std::cend(mapped_positions),
                                      original_mapping_position);
        }
    } // namespace
    
    double HaplotypeLikelihoodModel::log_probability(const AlignedRead& read, const Haplotype& haplotype,
                                                     const InactiveRegionState flank_state) const
    {
        Model model {2, 3}; // TODO: make part of an error model
        
        if (flank_state == InactiveRegionState::Unclear) model.flank_clear = false;
        
        auto gap_open_penalities = indel_error_model_.calculate_gap_open_penalties(haplotype);
        
        const auto best_mapping_positions = mapper_.map(read, haplotype);
        
        auto max_log_probability = std::numeric_limits<double>::lowest();
        
        for (const auto position : best_mapping_positions) {
            auto cur = compute_log_conditional_probability(haplotype.get_sequence(), read.get_sequence(),
                                                           read.get_qualities(), gap_open_penalities,
                                                           position, model);
            if (cur > max_log_probability) max_log_probability = cur;
        }
        
        assert(contains(haplotype, read));
        
        const auto original_mapping_position = begin_distance(read, haplotype);
        
        if (!contains(best_mapping_positions, original_mapping_position) && max_log_probability < 0) {
            auto cur = compute_log_conditional_probability(haplotype.get_sequence(), read.get_sequence(),
                                                           read.get_qualities(), gap_open_penalities,
                                                           original_mapping_position, model);
            if (cur > max_log_probability) max_log_probability = cur;
        }
        
        return max_log_probability;
    }
} // namespace Octopus
