//
//  haplotype_liklihood_model.cpp
//  Octopus
//
//  Created by Daniel Cooke on 25/09/2015.
//  Copyright © 2015 Oxford University. All rights reserved.
//

#include "haplotype_liklihood_model.hpp"

#include <utility>
#include <memory>
#include <cmath>
#include <limits>
#include <stdexcept>
#include <cassert>

#include "mappable.hpp"
#include "aligned_read.hpp"

#include <iostream> // TEST

namespace Octopus
{
void HaplotypeLikelihoodModel::set(const Haplotype& haplotype, boost::optional<FlankState> flank_state)
{
    haplotype_ = std::addressof(haplotype);
    haplotype_flank_state_ = std::move(flank_state);
    
    indel_error_model_.fill_gap_open_penalties(haplotype, haplotype_gap_open_penalities_);
    haplotype_gap_extension_penalty_ = indel_error_model_.calculate_gap_extension_penalty(haplotype);
}

void HaplotypeLikelihoodModel::clear() noexcept
{
    haplotype_ = nullptr;
    haplotype_flank_state_ = boost::none;
}

namespace
{
    std::size_t num_out_of_range_bases(const std::size_t mapping_position,
                                       const AlignedRead& read, const Haplotype& haplotype)
    {
        const auto alignment_size = sequence_size(read) + mapping_position + PairHMM::AlignmenetPad;
        
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

template <typename InputIt, typename T>
double log_probability(const AlignedRead& read, const Haplotype& haplotype,
                       InputIt first_mapping_position, InputIt last_mapping_position,
                       const T& gap_open_penalities, const PairHMM::Model& model)
{
    assert(contains(haplotype, read));
    
    const auto original_mapping_position = begin_distance(read, haplotype);
    
    auto max_log_probability = std::numeric_limits<double>::lowest();
    
    bool is_original_position_mapped {false}, has_in_range_mapping_position {false};
    
    std::for_each(first_mapping_position, last_mapping_position,
                  [&] (const auto position) {
                      if (position == original_mapping_position) {
                          is_original_position_mapped = true;
                      }
                      
                      if (is_in_range(position, read, haplotype)) {
                          has_in_range_mapping_position = true;
                          
                          auto cur = PairHMM::align_around_offset(haplotype.sequence(),
                                                                  read.sequence(),
                                                                  read.qualities(),
                                                                  gap_open_penalities,
                                                                  position,
                                                                  model);
                          
                          if (cur > max_log_probability) {
                              max_log_probability = cur;
                          }
                      }
                  });
    
    if (!is_original_position_mapped && is_in_range(original_mapping_position, read, haplotype)) {
        has_in_range_mapping_position = true;
        
        auto cur = PairHMM::align_around_offset(haplotype.sequence(), read.sequence(),
                                                read.qualities(), gap_open_penalities,
                                                original_mapping_position, model);
        
        if (cur > max_log_probability) {
            max_log_probability = cur;
        }
    }
    
    if (!has_in_range_mapping_position) {
        const auto min_shift = num_out_of_range_bases(original_mapping_position, read, haplotype);
        
        assert(original_mapping_position >= min_shift);
        
        const auto final_mapping_position = original_mapping_position - min_shift;
        
        max_log_probability = PairHMM::align_around_offset(haplotype.sequence(),
                                                           read.sequence(),
                                                           read.qualities(),
                                                           gap_open_penalities,
                                                           final_mapping_position,
                                                           model);
    }
    
    assert(max_log_probability > std::numeric_limits<double>::lowest());
    
    return max_log_probability;
}

// public methods

HaplotypeLikelihoodModel::HaplotypeLikelihoodModel()
: HaplotypeLikelihoodModel {2, ReadIndelErrorModel {}} {}

HaplotypeLikelihoodModel::HaplotypeLikelihoodModel(PenaltyType base_change_penalty,
                                                   ReadIndelErrorModel indel_model)
:
base_change_penalty_ {base_change_penalty},
indel_error_model_ {std::move(indel_model)},
haplotype_ {nullptr},
haplotype_flank_state_ {},
haplotype_gap_open_penalities_ {},
haplotype_gap_extension_penalty_ {}
{}

HaplotypeLikelihoodModel::HaplotypeLikelihoodModel(PenaltyType base_change_penalty,
                                                   ReadIndelErrorModel indel_model,
                                                   const Haplotype& haplotype,
                                                   boost::optional<FlankState> flank_state)
: HaplotypeLikelihoodModel {base_change_penalty, std::move(indel_model)}
{
    this->set(haplotype, std::move(flank_state));
}

double HaplotypeLikelihoodModel::log_probability(const AlignedRead& read,
                                                 MapPositionItr first_mapping_position,
                                                 MapPositionItr last_mapping_position) const
{
    if (haplotype_ == nullptr) {
        throw std::runtime_error {"HaplotypeLikelihoodModel: no buffered Haplotype"};
    }
    
    PairHMM::Model hmm_model;
    hmm_model.nucprior  = base_change_penalty_;
    hmm_model.gapextend = haplotype_gap_extension_penalty_;
    
    if (haplotype_flank_state_) {
        hmm_model.lhs_flank_size = haplotype_flank_state_->lhs_flank;
        hmm_model.rhs_flank_size = haplotype_flank_state_->rhs_flank;
    } else {
        hmm_model.lhs_flank_size = 0;
        hmm_model.rhs_flank_size = 0;
    }
    
    return Octopus::log_probability(read, *haplotype_,
                                    first_mapping_position, last_mapping_position,
                                    haplotype_gap_open_penalities_, hmm_model);
}
} // namespace Octopus
