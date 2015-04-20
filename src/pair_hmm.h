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

#include "aligned_read.h"
#include "maths.h"

using std::size_t;
using std::uint8_t;

struct RandomModel
{
    RandomModel() = default;
    
    double background_probability;
    double end_probability;
};

struct MatchModel
{
    MatchModel() = default;
    
    double match_probability;
    double gap_open_probability;
    double gap_extend_probability;
    double end_probability;
};

template <typename T>
struct LocalPairHmmLogState
{
    static constexpr const T Negative_infinity = std::numeric_limits<T>::lowest();
    
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
    
    LocalPairHmmLogState(T x)
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
    
    T begin;
    T random_x_1;
    T silent_1;
    T random_y_1;
    T silent_2;
    T match;
    T insertion;
    T deletion;
    T silent_3;
    T random_x_2;
    T silent_4;
    T random_y_2;
};

template <typename T>
T nuc_log_viterbi_local(const AlignedRead::SequenceType& sequence1,
                        const AlignedRead::SequenceType& sequence2,
                        const AlignedRead::Qualities& quals, MatchModel m, RandomModel r)
{
    static T log_prob_background = std::log(0.25);
    
    T log_prob_gap_open       = std::log(m.gap_open_probability);
    T log_prob_gap_extend     = std::log(m.gap_extend_probability);
    T log_prob_match_end      = std::log(m.end_probability);
    T log_prob_rand_end       = std::log(r.end_probability);
    T log_prob_rand_cont      = std::log(1 - r.end_probability);
    T log_prob_continue_match = std::log(1 - 2 * m.gap_open_probability - m.end_probability);
    T log_prob_to_match       = std::log(1 - m.gap_extend_probability - m.end_probability);
    
    static T ln_10_div_10 = std::log(10) / 10;
    static T ln_3         = std::log(3);
    
    auto sequence1_length = sequence1.size();
    auto sequence2_length = sequence2.size();
    
    std::vector<LocalPairHmmLogState<T>> current_column(sequence2_length + 2);
    std::vector<LocalPairHmmLogState<T>> previous_column(sequence2_length + 2);
    
    current_column[1].begin = 0;
    
    for (size_t i {1}; i <= sequence1_length + 1; ++i) {
        for (size_t j {1}; j <= sequence2_length + 1; ++j) {
            
            current_column[j].random_x_1 = log_prob_background + log_prob_rand_cont + std::max({
                previous_column[j].begin,
                previous_column[j].random_x_1
            });
            
            current_column[j].silent_1 = log_prob_rand_end + std::max({
                current_column[j].begin,
                current_column[j].random_x_1
            });
            
            current_column[j].random_y_1 = log_prob_background + log_prob_rand_cont + std::max({
                current_column[j - 1].silent_1,
                current_column[j - 1].random_y_1
            });
            
            current_column[j].silent_2 = log_prob_rand_end + std::max({
                current_column[j].silent_1,
                current_column[j].random_y_1
            });
            
            if (i > 1 && j > 1) {
                auto match_log_prob = ((sequence1[i - 2] == sequence2[j - 2]) ?
                                       std::log(1 - std::pow(10, -static_cast<T>(quals[j - 2]) / 10)) :
                                       -ln_10_div_10 * quals[j - 2] - ln_3);
                current_column[j].match = match_log_prob + std::max({
                    log_prob_continue_match + previous_column[j - 1].match,
                    log_prob_to_match + previous_column[j - 1].insertion,
                    log_prob_to_match + previous_column[j - 1].deletion,
                    log_prob_continue_match + previous_column[j - 1].silent_2
                });
            }
            
            current_column[j].insertion = log_prob_background + std::max({
                log_prob_gap_open + previous_column[j].match,
                log_prob_gap_extend + previous_column[j].insertion,
                log_prob_gap_open + previous_column[j].silent_2
            });
            
            current_column[j].deletion = log_prob_background + std::max({
                log_prob_gap_open + current_column[j - 1].match,
                log_prob_gap_extend + current_column[j - 1].insertion,
                log_prob_gap_open + current_column[j - 1].silent_2
            });
            
            current_column[j].silent_3 = log_prob_match_end + std::max({
                current_column[j].match,
                current_column[j].insertion,
                current_column[j].deletion,
                current_column[j].silent_2
            });
            
            current_column[j].random_x_2 = log_prob_background + log_prob_rand_cont + std::max({
                previous_column[j].silent_3,
                previous_column[j].random_x_2
            });
            
            current_column[j].silent_4 = log_prob_rand_end + std::max({
                current_column[j].silent_3,
                current_column[j].random_x_2
            });
            
            current_column[j].random_y_2 = log_prob_background + log_prob_rand_cont + std::max({
                current_column[j - 1].silent_4,
                current_column[j - 1].random_y_2
            });
        }
        
        if (i == 2) previous_column[1].begin = LocalPairHmmLogState<T>::Negative_infinity;
        
        std::swap(current_column, previous_column);
    }
    
    return log_prob_rand_end + std::max({
        previous_column[sequence2_length + 1].silent_4,
        previous_column[sequence2_length + 1].random_y_2
    });
}

template <typename T>
T nuc_log_forward_local(const std::string& sequence1, const std::string& sequence2,
                        const std::vector<uint8_t>& quals, MatchModel m, RandomModel r)
{
    static T log_prob_background = std::log(0.25);
    
    T log_prob_gap_open       = std::log(m.gap_open_probability);
    T log_prob_gap_extend     = std::log(m.gap_extend_probability);
    T log_prob_match_end      = std::log(m.end_probability);
    T log_prob_rand_end       = std::log(r.end_probability);
    T log_prob_rand_cont      = std::log(1 - r.end_probability);
    T log_prob_continue_match = std::log(1 - 2 * m.gap_open_probability - m.end_probability);
    T log_prob_to_match       = std::log(1 - m.gap_extend_probability - m.end_probability);
    
    static T ln_10_div_10 = std::log(10) / 10;
    static T ln_3         = std::log(3);
    
    auto sequence1_length = sequence1.size();
    auto sequence2_length = sequence2.size();
    
    std::vector<LocalPairHmmLogState<T>> current_column(sequence2_length + 2);
    std::vector<LocalPairHmmLogState<T>> previous_column(sequence2_length + 2);
    
    current_column[1].begin = 0;
    
    for (size_t i {1}; i <= sequence1_length + 1; ++i) {
        for (size_t j {1}; j <= sequence2_length + 1; ++j) {
            current_column[j].random_x_1 = log_prob_background + log_prob_rand_cont + log_sum_exp(
                previous_column[j].begin,
                previous_column[j].random_x_1
            );
            
            current_column[j].silent_1 = log_prob_rand_end + log_sum_exp(
                current_column[j].begin,
                current_column[j].random_x_1
            );
            
            current_column[j].random_y_1 = log_prob_background + log_prob_rand_cont + log_sum_exp(
                current_column[j - 1].silent_1,
                current_column[j - 1].random_y_1
            );
            
            current_column[j].silent_2 = log_prob_rand_end + log_sum_exp(
                current_column[j].silent_1,
                current_column[j].random_y_1
            );
            
            if (i > 1 && j > 1) {
                auto match_log_prob = ((sequence1[i - 2] == sequence2[j - 2]) ?
                                       std::log(1 - std::pow(10, -static_cast<T>(quals[j - 2]) / 10)) :
                                       -ln_10_div_10 * quals[j - 2] - ln_3);
                current_column[j].match = match_log_prob + log_sum_exp(
                    log_prob_continue_match + previous_column[j - 1].match,
                    log_prob_to_match + previous_column[j - 1].insertion,
                    log_prob_to_match + previous_column[j - 1].deletion,
                    log_prob_continue_match + previous_column[j - 1].silent_2
                );
            }
            
            current_column[j].insertion = log_prob_background + log_sum_exp(
                log_prob_gap_open + previous_column[j].match,
                log_prob_gap_extend + previous_column[j].insertion,
                log_prob_gap_open + previous_column[j].silent_2
            );
            
            current_column[j].deletion = log_prob_background + log_sum_exp(
                log_prob_gap_open + current_column[j - 1].match,
                log_prob_gap_extend + current_column[j - 1].insertion,
                log_prob_gap_open + current_column[j - 1].silent_2
            );
            
            current_column[j].silent_3 = log_prob_match_end + log_sum_exp(
                current_column[j].match,
                current_column[j].insertion,
                current_column[j].deletion,
                current_column[j].silent_2
            );
            
            current_column[j].random_x_2 = log_prob_background + log_prob_rand_cont + log_sum_exp(
                previous_column[j].silent_3,
                previous_column[j].random_x_2
            );
            
            current_column[j].silent_4 = log_prob_rand_end + log_sum_exp(
                current_column[j].silent_3,
                current_column[j].random_x_2
            );
            
            current_column[j].random_y_2 = log_prob_background + log_prob_rand_cont + log_sum_exp(
                current_column[j - 1].silent_4,
                current_column[j - 1].random_y_2
            );
        }
        
        if (i == 2) previous_column[1].begin = LocalPairHmmLogState<T>::Negative_infinity;
            
        std::swap(current_column, previous_column);
    }
            
    return log_prob_rand_end + log_sum_exp(
        previous_column[sequence2_length + 1].silent_4,
        previous_column[sequence2_length + 1].random_y_2
    );
}

#endif
