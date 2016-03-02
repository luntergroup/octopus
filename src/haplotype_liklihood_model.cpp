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
#include <cassert>
#include <limits>

#include "mappable.hpp"
#include "aligned_read.hpp"

#include <iostream> // TEST
#include <chrono>   // TEST

namespace Octopus
{
    // public methods
    
    HaplotypeLikelihoodModel::HaplotypeLikelihoodModel(const Haplotype& haplotype,
                                                       InactiveRegionState flank_state)
    :
    indel_error_model_ {},
    haplotype_ {haplotype},
    haplotype_gap_open_penalities_ {indel_error_model_.calculate_gap_open_penalties(haplotype)},
    model_ {2, 3, flank_state != InactiveRegionState::Clear}
    {}
    
    namespace
    {
        std::size_t num_out_of_range_bases(const std::size_t mapping_position,
                                           const AlignedRead& read, const Haplotype& haplotype)
        {
            const auto alignment_size = sequence_size(read) + mapping_position + 15;
            
            if (alignment_size > sequence_size(haplotype)) {
                return alignment_size - sequence_size(haplotype);
            }
            
            return 0;
        }
        
        bool is_in_range(const std::size_t mapping_position,
                             const AlignedRead& read, const Haplotype& haplotype)
        {
            return num_out_of_range_bases(mapping_position, read, haplotype) == 0;
        }
    } // namespace
    
    double HaplotypeLikelihoodModel::log_probability(const AlignedRead& read,
                                                     MapPositionItr first_mapping_position,
                                                     MapPositionItr last_mapping_position) const
    {
        return log_probability(read, haplotype_, first_mapping_position, last_mapping_position,
                               haplotype_gap_open_penalities_);
    }
    
    // private methods
    
    double HaplotypeLikelihoodModel::log_probability(const AlignedRead& read, const Haplotype& haplotype,
                                                     MapPositionItr first_mapping_position,
                                                     MapPositionItr last_mapping_position,
                                                     const std::vector<char>& gap_open_penalities) const
    {
        const auto original_mapping_position = begin_distance(read, haplotype);
        
        auto max_log_probability = std::numeric_limits<double>::lowest();
        
        bool is_original_position_mapped {false};
        bool has_in_range_mapping_position {false};
        
        std::for_each(first_mapping_position, last_mapping_position,
                      [&] (const auto position) {
                          if (is_in_range(position, read, haplotype)) {
                              has_in_range_mapping_position = true;
                              
                              auto cur = PairHMM::align_around_offset(haplotype.get_sequence(), read.get_sequence(),
                                                                      read.get_qualities(), gap_open_penalities,
                                                                      position, model_);
                              
                              if (cur > max_log_probability) {
                                  max_log_probability = cur;
                              }
                          }
                          
                          if (position == original_mapping_position) {
                              is_original_position_mapped = true;
                          }
                      });
        
        assert(contains(haplotype, read));
        
        if (!has_in_range_mapping_position) {
            auto final_mapping_position = original_mapping_position;
            
            if (is_original_position_mapped
                || !is_in_range(original_mapping_position, read, haplotype)) {
                const auto min_shift = num_out_of_range_bases(original_mapping_position, read, haplotype);
                assert(original_mapping_position >= min_shift);
                final_mapping_position -= min_shift;
            }
            
            max_log_probability = PairHMM::align_around_offset(haplotype.get_sequence(), read.get_sequence(),
                                                               read.get_qualities(), gap_open_penalities,
                                                               final_mapping_position, model_);
        }
        
        assert(max_log_probability > std::numeric_limits<double>::lowest());
        
        return max_log_probability;
    }
} // namespace Octopus
