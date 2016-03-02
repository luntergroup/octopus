//
//  pair_hmm.cpp
//  pair_hmm
//
//  Created by Daniel Cooke on 14/12/2015.
//  Copyright © 2015 Oxford University. All rights reserved.
//

#include "pair_hmm.hpp"

#include <stdexcept>
#include <algorithm>
#include <numeric>
#include <iterator>
#include <cmath>
#include <functional>
#include <type_traits>
#include <limits>
#include <array>
#include <cassert>

#include "align.h"
#include "banded_simd_viterbi.hpp"

#include <iostream> // DEBUG

namespace PairHMM
{
static constexpr double ln_10_div_10 {0.23025850929940458};

template <typename T, typename = std::enable_if_t<std::is_integral<T>::value>>
constexpr auto num_values() noexcept
{
    return std::numeric_limits<T>::max() - std::numeric_limits<T>::min() + 1;
}

template <std::size_t N>
auto make_phred_to_ln_prob_lookup() noexcept
{
    std::array<double, N> result {};
    
    for (std::size_t i {0}; i < N; ++i) {
        result[i] = -ln_10_div_10 * i;
    }
    
    return result;
}

template <typename InputIt1, typename InputIt2>
auto count_mismatches(InputIt1 first1, InputIt1 last1, InputIt2 first2)
{
    return std::inner_product(first1, last1, first2, std::size_t {0},
                              std::plus<void>(), std::not_equal_to<void>());
}

auto align(const std::string& truth, const std::string& target,
           const std::vector<std::uint8_t>& target_qualities,
           const std::vector<char>& truth_gap_open_penalties,
           const std::size_t target_offset, const Model& model)
{
    const auto truth_alignment_size = static_cast<int>(target.size() + 15);
    
    if (target_offset + truth_alignment_size > truth.size()) {
        return std::numeric_limits<double>::lowest();
    }
    
    if (!model.do_backtrace) {
        const auto score = fastAlignmentRoutine(truth.data() + target_offset,
                                                target.data(),
                                                reinterpret_cast<const char*>(target_qualities.data()),
                                                //truncated_target_qualities.data(),
                                                truth_alignment_size,
                                                static_cast<int>(target.size()),
                                                static_cast<int>(model.gapextend),
                                                static_cast<int>(model.nucprior),
                                                truth_gap_open_penalties.data() + target_offset);
        
        return -ln_10_div_10 * static_cast<double>(score);
    }
    
    int first_pos;
    std::vector<char> align1(2 * target.size() + 16), align2(2 * target.size() + 16);
    
    const auto score = fastAlignmentRoutine(truth.data() + target_offset,
                                            target.data(),
                                            reinterpret_cast<const char*>(target_qualities.data()),
                                            truth_alignment_size,
                                            static_cast<int>(target.size()),
                                            static_cast<int>(model.gapextend),
                                            static_cast<int>(model.nucprior),
                                            truth_gap_open_penalties.data() + target_offset,
                                            align1.data(), align2.data(), &first_pos);
    
    const auto flank_score = calculateFlankScore(truth_alignment_size,
                                                 0,
                                                 reinterpret_cast<const char*>(target_qualities.data()),
                                                 truth_gap_open_penalties.data(),
                                                 static_cast<int>(model.gapextend),
                                                 static_cast<int>(model.nucprior),
                                                 static_cast<int>(first_pos + target_offset),
                                                 align1.data(), align2.data());
    
    return -ln_10_div_10 * static_cast<double>(score - flank_score);
}

double align_around_offset(const std::string& truth, const std::string& target,
                           const std::vector<std::uint8_t>& target_qualities,
                           const std::vector<char>& truth_gap_open_penalties,
                           const std::size_t target_offset, const Model& model)
{
    using std::cbegin; using std::cend; using std::next;
    
    static constexpr std::size_t NUM_QUALITIES {num_values<std::uint8_t>()};
    
    static const auto phred_to_ln_probability = make_phred_to_ln_prob_lookup<NUM_QUALITIES>();
    
    assert(target.size() == target_qualities.size());
    assert(truth.size() == truth_gap_open_penalties.size());
    assert(std::max(truth.size(), target.size()) > target_offset);
    
    if (target_offset + target.size() > truth.size()) {
        return std::numeric_limits<double>::lowest();
    }
    
    const auto offsetted_truth_begin_itr = next(cbegin(truth), target_offset);
    
    const auto p = std::mismatch(cbegin(target), cend(target), offsetted_truth_begin_itr);
    
    if (p.first == cend(target)) return 0;
    
    const auto num_mismatches = count_mismatches(next(p.first), cend(target), next(p.second)) + 1;
    
    if (num_mismatches == 1) {
        const auto mismatch_index = std::distance(offsetted_truth_begin_itr, p.second);
        
        if (target_qualities[mismatch_index] <= truth_gap_open_penalties[mismatch_index]) {
            return phred_to_ln_probability[target_qualities[mismatch_index]];
        }
        
        if (std::equal(next(p.first), cend(target), p.second)) {
            return phred_to_ln_probability[truth_gap_open_penalties[mismatch_index]];
        }
        
        return phred_to_ln_probability[target_qualities[mismatch_index]];
    }
    
    return align(truth, target, target_qualities, truth_gap_open_penalties,
                 target_offset, model);
    
    return 0;
}
    
} // namespace PairHMM
