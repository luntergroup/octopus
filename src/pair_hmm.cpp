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

#include "align.h"
#include "banded_simd_viterbi.hpp"

#include <iostream> // DEBUG
using std::cout;
using std::endl;

namespace
{
    static constexpr double ln_10_div_10 {0.23025850929940458};
    
    template <typename T, typename = std::enable_if_t<std::is_integral<T>::value>>
    constexpr std::size_t num_values() noexcept
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
               const std::vector<std::uint8_t>& target_gap_open_penalties,
               const std::size_t offset_hint,
               const Model& model)
    {
        const auto score = fastAlignmentRoutine(truth.c_str() + offset_hint,
                                                target.c_str(),
                                                reinterpret_cast<const char*>(target_qualities.data()),
                                                static_cast<int>(target.size() + 15),
                                                static_cast<int>(target.size()),
                                                static_cast<int>(model.gapextend),
                                                static_cast<int>(model.nucprior),
                                                reinterpret_cast<const char*>(target_gap_open_penalties.data()) + offset_hint,
                                                nullptr, nullptr, nullptr);
        
        //    const auto score = align(truth.c_str() + offset_hint, target.c_str(),
        //                             reinterpret_cast<const char*>(target_qualities.data()),
        //                             static_cast<unsigned>(truth.size()),
        //                             static_cast<unsigned>(target.size()),
        //                             static_cast<std::uint8_t>(model.gapextend),
        //                             static_cast<std::uint8_t>(model.nucprior),
        //                             reinterpret_cast<const std::uint8_t*>(target_gap_open_penalties.data()));
        
        return -ln_10_div_10 * static_cast<double>(score);
    }
} // namespace

double compute_log_conditional_probability(const std::string& truth, const std::string& target,
                                           const std::vector<std::uint8_t>& target_qualities,
                                           const std::vector<std::uint8_t>& target_gap_open_penalties,
                                           const std::size_t target_offset_into_truth_hint,
                                           const Model& model)
{
    using std::cbegin; using std::cend; using std::next;
    
    static constexpr std::size_t NUM_QUALITIES {num_values<std::uint8_t>()};
    
    static const auto phred_to_ln_probability = make_phred_to_ln_prob_lookup<NUM_QUALITIES>();
    
    if (target.size() != target_qualities.size()) {
        throw std::logic_error {"compute_log_conditional_probability: target size does not match"
                                " target_qualities size"};
    }
    
    if (truth.size() != target_gap_open_penalties.size()) {
        throw std::logic_error {"compute_log_conditional_probability: truth size does not match"
                                " target_gap_open_penalties size"};
    }
    
    if (std::max(truth.size(), target.size()) <= target_offset_into_truth_hint) {
        throw std::logic_error {"compute_log_conditional_probability: offset hint is not"
                                " in range"};
    }
    
    const auto hinted_truth_begin_itr = next(cbegin(truth), target_offset_into_truth_hint);
    
    const auto p = std::mismatch(cbegin(target), cend(target), hinted_truth_begin_itr);
    
    if (p.first == cend(target)) {
        return 0;
    }
    
    const auto num_mismatches = count_mismatches(next(p.first), cend(target), next(p.second)) + 1;
    
    if (num_mismatches == 1) {
        const auto mismatch_index = std::distance(hinted_truth_begin_itr, p.second);
        
        if (target_qualities[mismatch_index] <= target_gap_open_penalties[mismatch_index]) {
            return phred_to_ln_probability[target_qualities[mismatch_index]];
        }
        
        if (std::equal(next(p.first), cend(target), p.second)) {
            return phred_to_ln_probability[target_gap_open_penalties[mismatch_index]];
        }
        
        return phred_to_ln_probability[target_qualities[mismatch_index]];
    }
    
    return align(truth, target, target_qualities, target_gap_open_penalties,
                 target_offset_into_truth_hint, model);
    
    return 0;
}
