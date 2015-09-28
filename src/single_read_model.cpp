//
//  single_read_model.cpp
//  Octopus
//
//  Created by Daniel Cooke on 25/09/2015.
//  Copyright © 2015 Oxford University. All rights reserved.
//

#include "single_read_model.hpp"

#include <cmath>     // std::log
#include <algorithm> // std::max, std::min

#include "mappable.hpp"
#include "pair_hmm.hpp"

#include <iostream> // TEST

namespace Octopus
{
    // public methods
    
    SingleReadModel::SingleReadModel(size_t max_cache_size)
    :
    cache_ {},
    max_cache_size_ {max_cache_size},
    cache_size_ {}
    {
        cache_.reserve(max_cache_size / 2);
    }
    
    SingleReadModel::SingleReadModel(size_t max_num_reads, size_t max_num_haplotypes)
    :
    cache_ {},
    max_cache_size_ {max_num_reads * max_num_haplotypes},
    cache_size_ {}
    {
        cache_.reserve(max_num_reads);
    }
    
    SingleReadModel::RealType SingleReadModel::log_probability(const AlignedRead& read, const Haplotype& haplotype)
    {
        if (is_cached(read, haplotype)) return get_cached(read, haplotype);
        
        // TODO: make these members when pair_hmm is finalised
        
        RandomModel<RealType> lhs_random {};
        lhs_random.target_emission_probability = 0.25;
        lhs_random.query_emission_probability  = 0.25;
        
        MatchModel<RealType> match {};
        match.match_probability      = 0.25;
        match.gap_open_probability   = 0.015; // TODO: should be part of an error model
        match.gap_extend_probability = 0.020; // TODO: should be part of an error model
        
        RandomModel<RealType> rhs_random {};
        rhs_random.target_emission_probability = 0.25;
        rhs_random.query_emission_probability  = 0.25;
        
        // m.end_probability must satisfy:
        // m.end_probability <= 1 - 2 * m.gap_open_probability
        // m.end_probability <= 1 - m.gap_extend_probability
        
        auto max_match_end_prob = 1 - std::max(2 * match.gap_open_probability, match.gap_extend_probability);
        
        if (overlaps(read, haplotype)) {
            auto overlapped_region = get_overlapped(read, haplotype);
            auto covered_region    = get_encompassing(read, haplotype);
            
            if (begins_before(read, haplotype)) {
                lhs_random.target_end_probability = 0.99;
                lhs_random.query_end_probability  = 1.0 / (size(get_left_overhang(covered_region, overlapped_region)) + 1);
            } else {
                lhs_random.target_end_probability = 1.0 / (size(get_left_overhang(covered_region, overlapped_region)) + 1);
                lhs_random.query_end_probability  = 0.99;
            }
            
            match.end_probability = std::min(1.0 / (size(overlapped_region) + 1), max_match_end_prob);
            
            if (ends_before(read, haplotype)) {
                rhs_random.target_end_probability = 1.0 / (size(get_right_overhang(covered_region, overlapped_region)) + 1);
                rhs_random.query_end_probability  = 0.99;
            } else {
                rhs_random.target_end_probability = 0.99;
                rhs_random.query_end_probability  = 1.0 / (size(get_right_overhang(covered_region, overlapped_region)) + 1);
            }
        } else {
            lhs_random.target_end_probability = 1.0 / (size(haplotype) + 1);
            lhs_random.query_end_probability  = 1.0 / (size(read) + 1);
            
            match.end_probability = max_match_end_prob;
            
            rhs_random.target_end_probability = 0.99;
            rhs_random.query_end_probability  = 0.99;
        }
        
        auto joint_log_probability = nuc_log_viterbi_local<RealType>(haplotype.get_sequence(), read.get_sequence(),
                                                                     read.get_qualities(),
                                                                     lhs_random, match, rhs_random);
        
        auto conditional_log_probability = joint_log_probability -
        haplotype.get_sequence().size() * std::log(lhs_random.target_emission_probability);
        
        cache(read, haplotype, conditional_log_probability);
        
        return conditional_log_probability;
    }
    
    void SingleReadModel::clear_cache()
    {
        cache_.clear();
        cache_size_ = 0;
    }
    
    // private methods
    
    bool SingleReadModel::is_cached(const AlignedRead& read, const Haplotype& haplotype) const noexcept
    {
        return cache_.count(read) == 1 && cache_.at(read).count(haplotype) == 1;
    }
    
    void SingleReadModel::cache(const AlignedRead& read, const Haplotype& haplotype, RealType value)
    {
        if (cache_size_ < max_cache_size_) {
            cache_[read][haplotype] = value;
            ++cache_size_;
        }
    }
    
    SingleReadModel::RealType SingleReadModel::get_cached(const AlignedRead& read, const Haplotype& haplotype) const
    {
        return cache_.at(read).at(haplotype);
    }

} // namespace Octopus
