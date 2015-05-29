//
//  pair_hmm.h
//  Octopus
//
//  Created by Daniel Cooke on 27/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_pair_hmm_h
#define Octopus_pair_hmm_h

#include <vector>
#include <string>
#include <cstddef>
#include <algorithm>
#include <cmath>
#include <limits>
#include <cstdint>
#include <limits>
#include <numeric>
#include <array>

#include "maths.h"

using std::size_t;

template <typename RealType>
struct RandomModel
{
    RandomModel() = default;
    
    RealType target_emission_probability;
    RealType query_emission_probability;
    RealType target_end_probability;
    RealType query_end_probability;
};

template <typename RealType>
struct MatchModel
{
    MatchModel() = default;
    
    RealType match_probability;
    RealType gap_open_probability;
    RealType gap_extend_probability;
    RealType end_probability;
};

template <typename RealType>
struct LocalPairHmmLogState
{
    static constexpr const RealType Negative_infinity = std::numeric_limits<RealType>::lowest();
    
    LocalPairHmmLogState()
    :
    begin {Negative_infinity},
    random_x_1 {Negative_infinity},
    silent_1 {Negative_infinity},
    random_y_1 {Negative_infinity},
    silent_2 {Negative_infinity},
    match {Negative_infinity},
    insertion {Negative_infinity},
    deletion {Negative_infinity},
    silent_3 {Negative_infinity},
    random_x_2 {Negative_infinity},
    silent_4 {Negative_infinity},
    random_y_2 {Negative_infinity}
    {}
    
    LocalPairHmmLogState(RealType x)
    :
    begin {x},
    random_x_1 {x},
    silent_1 {x},
    random_y_1 {x},
    silent_2 {x},
    match {x},
    insertion {x},
    deletion {x},
    silent_3 {x},
    random_x_2 {x},
    silent_4 {x},
    random_y_2 {x}
    {}
    
    RealType begin;
    RealType random_x_1;
    RealType silent_1;
    RealType random_y_1;
    RealType silent_2;
    RealType match;
    RealType insertion;
    RealType deletion;
    RealType silent_3;
    RealType random_x_2;
    RealType silent_4;
    RealType random_y_2;
};

template <typename RealType>
std::vector<RealType> make_log_prob_match_lookup()
{
    static constexpr const unsigned Num_qualities {256};
    
    std::vector<RealType> the_lookup(Num_qualities);
    
    for (unsigned i {}; i < Num_qualities; ++i) {
        the_lookup[i] = std::log(1 - std::pow(10, -static_cast<RealType>(i) / 10));
    }
    
    return the_lookup;
}

template <typename RealType>
std::vector<RealType> make_log_prob_mismatch_lookup()
{
    static constexpr const unsigned Num_qualities {256};
    
    static const RealType ln_10_div_10 {static_cast<RealType>(std::log(10)) / 10};
    static const RealType ln_3         {static_cast<RealType>(std::log(3))};
    
    std::vector<RealType> the_lookup(Num_qualities);
    
    for (unsigned i {}; i < Num_qualities; ++i) {
        the_lookup[i] = -ln_10_div_10 * i - ln_3;
    }
    
    return the_lookup;
}

template <typename RealType, typename SequenceType1, typename SequenceType2>
RealType nuc_log_viterbi_local(const SequenceType1& target, const SequenceType1& query,
                               const SequenceType2& query_qualities, 
                               MatchModel<RealType> m, RandomModel<RealType> r1, 
                               RandomModel<RealType> r2)
{
    RealType log_prob_target_rand1      = std::log(r1.target_emission_probability);
    RealType log_prob_query_rand1       = std::log(r1.query_emission_probability);
    RealType log_prob_target_rand1_end  = std::log(r1.target_end_probability);
    RealType log_prob_target_rand1_cont = std::log(1 - r1.target_end_probability);
    RealType log_prob_query_rand1_end   = std::log(r1.target_end_probability);
    RealType log_prob_query_rand1_cont  = std::log(1 - r1.target_end_probability);
    
    RealType log_prob_match             = std::log(m.match_probability);
    RealType log_prob_gap_open          = std::log(m.gap_open_probability);
    RealType log_prob_gap_extend        = std::log(m.gap_extend_probability);
    RealType log_prob_match_end         = std::log(m.end_probability);
    RealType log_prob_continue_match    = std::log(1 - 2 * m.gap_open_probability - m.end_probability);
    RealType log_prob_to_match          = std::log(1 - m.gap_extend_probability - m.end_probability);
    
    RealType log_prob_target_rand2      = std::log(r2.target_emission_probability);
    RealType log_prob_query_rand2       = std::log(r2.query_emission_probability);
    RealType log_prob_target_rand2_end  = std::log(r2.target_end_probability);
    RealType log_prob_target_rand2_cont = std::log(1 - r2.target_end_probability);
    RealType log_prob_query_rand2_end   = std::log(r2.target_end_probability);
    RealType log_prob_query_rand2_cont  = std::log(1 - r2.target_end_probability);
    
    static auto log_prob_match_lookup    = make_log_prob_match_lookup<RealType>();
    static auto log_prob_mismatch_lookup = make_log_prob_mismatch_lookup<RealType>();
    
    auto target_length = target.size();
    auto query_length  = query.size();
    
    std::vector<LocalPairHmmLogState<RealType>> current_column(query_length + 2);
    std::vector<LocalPairHmmLogState<RealType>> previous_column(query_length + 2);
    
    current_column[1].begin = 0;
    
    for (size_t i {1}; i <= target_length + 1; ++i) {
        for (size_t j {1}; j <= query_length + 1; ++j) {
            
            current_column[j].random_x_1 = log_prob_target_rand1 + log_prob_target_rand1_cont + std::max(
                previous_column[j].begin,
                previous_column[j].random_x_1
            );
            
            current_column[j].silent_1 = log_prob_target_rand1_end + std::max(
                current_column[j].begin,
                current_column[j].random_x_1
            );
            
            current_column[j].random_y_1 = log_prob_query_rand1 + log_prob_query_rand1_cont + std::max(
                current_column[j - 1].silent_1,
                current_column[j - 1].random_y_1
            );
            
            current_column[j].silent_2 = log_prob_query_rand1_end + std::max(
                current_column[j].silent_1,
                current_column[j].random_y_1
            );
            
            if (i > 1 && j > 1) {
                current_column[j].match = log_prob_match + ((target[i - 2] == query[j - 2]) ?
                                           log_prob_match_lookup[query_qualities[j - 2]] :
                                           log_prob_mismatch_lookup[query_qualities[j - 2]])
                        + std::max({
                            log_prob_continue_match + previous_column[j - 1].match,
                            log_prob_to_match       + previous_column[j - 1].insertion,
                            log_prob_to_match       + previous_column[j - 1].deletion,
                            log_prob_continue_match + previous_column[j - 1].silent_2
                        });
            }
            
            current_column[j].insertion = log_prob_target_rand1 + std::max({
                log_prob_gap_open   + previous_column[j].match,
                log_prob_gap_extend + previous_column[j].insertion,
                log_prob_gap_open   + previous_column[j].silent_2
            });
            
            current_column[j].deletion = log_prob_query_rand1 + std::max({
                log_prob_gap_open   + current_column[j - 1].match,
                log_prob_gap_extend + current_column[j - 1].insertion,
                log_prob_gap_open   + current_column[j - 1].silent_2
            });
            
            current_column[j].silent_3 = log_prob_match_end + std::max({
                current_column[j].match,
                current_column[j].insertion,
                current_column[j].deletion,
                current_column[j].silent_2
            });
            
            current_column[j].random_x_2 = log_prob_target_rand2 + log_prob_target_rand2_cont + std::max(
                previous_column[j].silent_3,
                previous_column[j].random_x_2
            );
            
            current_column[j].silent_4 = log_prob_target_rand2_end + std::max(
                current_column[j].silent_3,
                current_column[j].random_x_2
            );
            
            current_column[j].random_y_2 = log_prob_query_rand2 + log_prob_query_rand2_cont + std::max(
                current_column[j - 1].silent_4,
                current_column[j - 1].random_y_2
            );
        }
        
        if (i == 2) previous_column[1].begin = LocalPairHmmLogState<RealType>::Negative_infinity;
        
        std::swap(current_column, previous_column);
    }
    
    return log_prob_query_rand2_end + std::max(
        previous_column[query_length + 1].silent_4,
        previous_column[query_length + 1].random_y_2
    );
}

//template <typename RealType, typename SequenceType1, typename SequenceType2>
//RealType nuc_log_forward_local(const SequenceType1& target, const SequenceType1& query,
//                               const SequenceType2& query_qualities, 
//                               MatchModel<RealType> m, RandomModel<RealType> r1, 
//                               RandomModel<RealType> r2)
//{
//    static auto log_prob_background = static_cast<RealType>(std::log(0.25));
//    
//    RealType log_prob_rand1_end      = std::log(r1.end_probability);
//    RealType log_prob_rand1_cont     = std::log(1 - r1.end_probability);
//    
//    RealType log_prob_gap_open       = std::log(m.gap_open_probability);
//    RealType log_prob_gap_extend     = std::log(m.gap_extend_probability);
//    RealType log_prob_match_end      = std::log(m.end_probability);
//    RealType log_prob_continue_match = std::log(1 - 2 * m.gap_open_probability - m.end_probability);
//    RealType log_prob_to_match       = std::log(1 - m.gap_extend_probability - m.end_probability);
//    
//    RealType log_prob_rand2_end      = std::log(r2.end_probability);
//    RealType log_prob_rand2_cont     = std::log(1 - r2.end_probability);
//    
//    static auto log_prob_match_lookup    = make_log_prob_match_lookup<RealType>();
//    static auto log_prob_mismatch_lookup = make_log_prob_mismatch_lookup<RealType>();
//    
//    auto target_length = target.size();
//    auto query_length  = query.size();
//    
//    std::vector<LocalPairHmmLogState<RealType>> current_column(query_length + 2);
//    std::vector<LocalPairHmmLogState<RealType>> previous_column(query_length + 2);
//    
//    current_column[1].begin = 0;
//    
//    for (size_t i {1}; i <= target_length + 1; ++i) {
//        for (size_t j {1}; j <= query_length + 1; ++j) {
//            
//            current_column[j].random_x_1 = log_prob_background + log_prob_rand1_cont + log_sum_exp(
//                previous_column[j].begin,
//                previous_column[j].random_x_1
//            );
//            
//            current_column[j].silent_1 = log_prob_rand1_end + log_sum_exp(
//                current_column[j].begin,
//                current_column[j].random_x_1
//            );
//            
//            current_column[j].random_y_1 = log_prob_background + log_prob_rand1_cont + log_sum_exp(
//                current_column[j - 1].silent_1,
//                current_column[j - 1].random_y_1
//            );
//            
//            current_column[j].silent_2 = log_prob_rand1_end + log_sum_exp(
//                current_column[j].silent_1,
//                current_column[j].random_y_1
//            );
//            
//            if (i > 1 && j > 1) {
//                current_column[j].match = ((target[i - 2] == query[j - 2]) ?
//                                           log_prob_match_lookup[query_qualities[j - 2]] :
//                                           log_prob_mismatch_lookup[query_qualities[j - 2]])
//                            + log_sum_exp(
//                                log_prob_continue_match + previous_column[j - 1].match,
//                                log_prob_to_match       + previous_column[j - 1].insertion,
//                                log_prob_to_match       + previous_column[j - 1].deletion,
//                                log_prob_continue_match + previous_column[j - 1].silent_2
//                            );
//            }
//            
//            current_column[j].insertion = log_prob_background + log_sum_exp(
//                log_prob_gap_open   + previous_column[j].match,
//                log_prob_gap_extend + previous_column[j].insertion,
//                log_prob_gap_open   + previous_column[j].silent_2
//            );
//            
//            current_column[j].deletion = log_prob_background + log_sum_exp(
//                log_prob_gap_open   + current_column[j - 1].match,
//                log_prob_gap_extend + current_column[j - 1].insertion,
//                log_prob_gap_open   + current_column[j - 1].silent_2
//            );
//            
//            current_column[j].silent_3 = log_prob_match_end + log_sum_exp(
//                current_column[j].match,
//                current_column[j].insertion,
//                current_column[j].deletion,
//                current_column[j].silent_2
//            );
//            
//            current_column[j].random_x_2 = log_prob_background + log_prob_rand2_cont + log_sum_exp(
//                previous_column[j].silent_3,
//                previous_column[j].random_x_2
//            );
//            
//            current_column[j].silent_4 = log_prob_rand2_end + log_sum_exp(
//                current_column[j].silent_3,
//                current_column[j].random_x_2
//            );
//            
//            current_column[j].random_y_2 = log_prob_background + log_prob_rand2_cont + log_sum_exp(
//                current_column[j - 1].silent_4,
//                current_column[j - 1].random_y_2
//            );
//        }
//        
//        if (i == 2) previous_column[1].begin = LocalPairHmmLogState<RealType>::Negative_infinity;
//          
//        std::swap(current_column, previous_column);
//    }
//            
//    return log_prob_rand2_end + log_sum_exp(
//        previous_column[query_length + 1].silent_4,
//        previous_column[query_length + 1].random_y_2
//    );
//}

#endif
