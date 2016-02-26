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
        std::size_t count_out_of_range_bases(const std::size_t mapping_position,
                                             const AlignedRead& read, const Haplotype& haplotype)
        {
            const auto alignment_size = sequence_size(read) + mapping_position + 15;
            
            if (alignment_size > sequence_size(haplotype)) {
                return alignment_size - sequence_size(haplotype);
            }
            
            return 0;
        }
        
        bool is_out_of_range(const std::size_t mapping_position,
                             const AlignedRead& read, const Haplotype& haplotype)
        {
            return count_out_of_range_bases(mapping_position, read, haplotype) > 0;
        }
        
        void remove_out_of_range_mapping_positions(std::vector<std::size_t>& mapping_positions,
                                                   const AlignedRead& read, const Haplotype& haplotype)
        {
            const auto it = std::remove_if(std::begin(mapping_positions), std::end(mapping_positions),
                                           [&] (const auto position) {
                                               return is_out_of_range(position, read, haplotype);
                                           });
            mapping_positions.erase(it, std::end(mapping_positions));
        }
        
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
        
        const auto gap_open_penalities = indel_error_model_.calculate_gap_open_penalties(haplotype);
        
        auto mapping_positions = mapper_.map(read, haplotype);
        
        remove_out_of_range_mapping_positions(mapping_positions, read, haplotype);
        
        const auto original_mapping_position = begin_distance(read, haplotype);
        
        if (mapping_positions.empty() && is_out_of_range(original_mapping_position, read, haplotype)) {
            const auto min_shift = count_out_of_range_bases(original_mapping_position, read, haplotype);
            assert(original_mapping_position >= min_shift);
            return compute_log_conditional_probability(haplotype.get_sequence(), read.get_sequence(),
                                                       read.get_qualities(), gap_open_penalities,
                                                       original_mapping_position - min_shift, model);
        }
        
        auto max_log_probability = std::numeric_limits<double>::lowest();
        
        for (const auto position : mapping_positions) {
            auto cur = compute_log_conditional_probability(haplotype.get_sequence(), read.get_sequence(),
                                                           read.get_qualities(), gap_open_penalities,
                                                           position, model);
            if (cur > max_log_probability) max_log_probability = cur;
        }
        
        assert(contains(haplotype, read));
        
        if (!contains(mapping_positions, original_mapping_position) && max_log_probability < 0) {
            auto cur = compute_log_conditional_probability(haplotype.get_sequence(), read.get_sequence(),
                                                           read.get_qualities(), gap_open_penalities,
                                                           original_mapping_position, model);
            if (cur > max_log_probability) max_log_probability = cur;
        }
        
        assert(max_log_probability > std::numeric_limits<double>::lowest());
        
        return max_log_probability;
    }
} // namespace Octopus
