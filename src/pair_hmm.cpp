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
#include <iostream>

#include "simd_pair_hmm.hpp"
#include "maths.hpp"
#include "cigar_string.hpp"

namespace octopus { namespace hmm {

using octopus::maths::constants::ln_10_div_10;

template <typename T, typename = std::enable_if_t<std::is_integral<T>::value>>
constexpr auto num_values() noexcept
{
    return std::numeric_limits<T>::max() - std::numeric_limits<T>::min() + 1;
}

template <std::size_t N>
constexpr auto make_phred_to_ln_prob_lookup() noexcept
{
    std::array<double, N> result {};
    
    for (std::size_t i {0}; i < N; ++i) {
        // This const_cast mess is because std::array::operator[] is not marked constexpr (until C++17)
        const_cast<double&>(static_cast<std::array<double, N> const&>(result)[i]) = -ln_10_div_10 * i;
    }
    
    return result;
}

template <typename T>
constexpr auto make_phred_to_ln_prob_lookup() noexcept
{
    return make_phred_to_ln_prob_lookup<num_values<T>()>();
}

bool is_target_in_truth_flank(const std::string& truth, const std::string& target,
                              const std::size_t target_offset, const Model& model)
{
    return target_offset < model.lhs_flank_size
                || (target_offset + target.size()) > (truth.size() - model.rhs_flank_size);
}

auto make_cigar(const std::vector<char>& align1, const std::vector<char>& align2)
{
    assert(!align1.empty());
    
    CigarString result {};
    result.reserve(align1.size());
    
    auto it1 = std::cbegin(align1);
    auto it2 = std::cbegin(align2);
    
    const auto last1 = std::find_if_not(std::crbegin(align1), std::crend(align1),
                                        [] (const auto x) { return x == 0; }).base();
    const auto last2 = std::next(it1, std::distance(it1, last1));
    
    while (it1 != last1) {
        const auto p = std::mismatch(it1, last1, it2);
        
        if (p.first != it1) {
            result.emplace_back(std::distance(it1, p.first), CigarOperation::SEQUENCE_MATCH);
            
            if (p.first == last1) break;
        }
        
        const static auto is_gap = [] (const auto b) { return b == '-'; };
        
        if (*p.first == '-') {
            const auto it3 = std::find_if_not(std::next(p.first), last1, is_gap);
            
            const auto n = std::distance(p.first, it3);
            
            result.emplace_back(n, CigarOperation::INSERTION);
            
            it1 = it3;
            it2 = std::next(p.second, n);
        } else if (*p.second == '-') {
            const auto it3 = std::find_if_not(std::next(p.second), last2, is_gap);
            
            const auto n = std::distance(p.second, it3);
            
            result.emplace_back(n, CigarOperation::DELETION);
            
            it1 = std::next(p.first, n);
            it2 = it3;
        } else {
            const auto p2 = std::mismatch(std::next(p.first), last1, std::next(p.second),
                                          std::not_equal_to<> {});
            
            result.emplace_back(std::distance(p.first, p2.first), CigarOperation::SUBSTITUTION);
            
            std::tie(it1, it2) = p2;
        }
    }
    
    result.shrink_to_fit();
    
    return result;
}

namespace debug {
    void print_alignment(const std::vector<char>& align1, const std::vector<char>& align2)
    {
        const auto isnt_null = [] (const char c) { return c != '\0'; };
        std::copy_if(std::cbegin(align1), std::cend(align1), std::ostreambuf_iterator<char>(std::cout),
                     isnt_null);
        std::cout << '\n';
        std::copy_if(std::cbegin(align2), std::cend(align2), std::ostreambuf_iterator<char>(std::cout),
                     isnt_null);
        std::cout << '\n';
    }
} // namespace debug

auto simd_align(const std::string& truth, const std::string& target,
                const std::vector<std::uint8_t>& target_qualities,
                const std::size_t target_offset,
                const Model& model)
{
    constexpr auto Pad = static_cast<int>(simd::min_flank_pad());
    
    const auto truth_alignment_size = static_cast<int>(target.size() + 2 * Pad - 1);
    
    const auto alignment_offset = static_cast<std::size_t>(std::max(0, static_cast<int>(target_offset) - Pad));
    
    if (alignment_offset + truth_alignment_size > truth.size()) {
        return std::numeric_limits<double>::lowest();
    }
    
    const auto qualities = reinterpret_cast<const std::int8_t*>(target_qualities.data());
    
    if (!is_target_in_truth_flank(truth, target, target_offset, model)) {
        const auto score = simd::align(truth.data() + alignment_offset,
                                       target.data(),
                                       qualities,
                                       truth_alignment_size, static_cast<int>(target.size()),
                                       model.snv_mask.data() + alignment_offset,
                                       model.snv_priors.data() + alignment_offset,
                                       model.gap_open_penalties.data() + alignment_offset,
                                       model.gap_extend, model.nuc_prior);
        
        return -ln_10_div_10 * static_cast<double>(score);
    }
    
    thread_local std::vector<char> align1 {}, align2 {};
    
    const auto max_alignment_size = 2 * (target.size() + Pad);
    
    align1.assign(max_alignment_size + 1, 0);
    align2.assign(max_alignment_size + 1, 0);
    
    int first_pos;
    
    const auto score = simd::align(truth.data() + alignment_offset,
                                   target.data(),
                                   qualities,
                                   truth_alignment_size,
                                   static_cast<int>(target.size()),
                                   model.snv_mask.data() + alignment_offset,
                                   model.snv_priors.data() + alignment_offset,
                                   model.gap_open_penalties.data() + alignment_offset,
                                   model.gap_extend, model.nuc_prior,
                                   align1.data(), align2.data(), first_pos);
    
//    debug::print_alignment(align1, align2);
//    std::cout << make_cigar(align1, align2) << std::endl;
    
    auto lhs_flank_size = static_cast<int>(model.lhs_flank_size);
    
    if (lhs_flank_size < alignment_offset) {
        lhs_flank_size = 0;
    } else {
        lhs_flank_size -= alignment_offset;
    }
    
    auto rhs_flank_size = static_cast<int>(model.rhs_flank_size);
    
    if (alignment_offset + truth_alignment_size < (truth.size() - model.rhs_flank_size)) {
        rhs_flank_size = 0;
    } else {
        rhs_flank_size += alignment_offset + truth_alignment_size;
        rhs_flank_size -= truth.size();
    }
    
    assert(lhs_flank_size >= 0 && rhs_flank_size >= 0);
    assert(align1.back() == 0); // required by calculate_flank_score
    
    const auto flank_score = simd::calculate_flank_score(truth_alignment_size,
                                                         lhs_flank_size, rhs_flank_size,
                                                         target.data(), qualities,
                                                         model.snv_mask.data() + alignment_offset,
                                                         model.snv_priors.data() + alignment_offset,
                                                         model.gap_open_penalties.data() + alignment_offset,
                                                         model.gap_extend, model.nuc_prior,
                                                         first_pos,
                                                         align1.data(), align2.data());
    
    assert(flank_score <= score);
    
    return -ln_10_div_10 * static_cast<double>(score - flank_score);
}

unsigned min_flank_pad() noexcept
{
    return simd::min_flank_pad();
}

void validate(const std::string& truth, const std::string& target,
              const std::vector<std::uint8_t>& target_qualities,
              const std::size_t target_offset,
              const Model& model)
{
    if (target.size() != target_qualities.size()) {
        throw std::invalid_argument {"PairHMM::align: target size not equal to target qualities length"};
    }
    if (truth.size() != model.snv_priors.size()) {
        throw std::invalid_argument {"PairHMM::align: truth size not equal to snv priors length"};
    }
    if (truth.size() != model.gap_open_penalties.size()) {
        throw std::invalid_argument {"PairHMM::align: truth size not equal to gap open penalties length"};
    }
    if (target_offset + target.size() > truth.size()) {
        throw std::invalid_argument {"PairHMM::align: target is not contained by truth"};
    }
}

double score(const std::string& truth, const std::string& target,
             const std::vector<std::uint8_t>& target_qualities,
             const std::size_t target_offset,
             const Model& model)
{
    using std::cbegin; using std::cend; using std::next; using std::distance;
    
    static constexpr auto Ln_probability = make_phred_to_ln_prob_lookup<std::uint8_t>();
    
    validate(truth, target, target_qualities, target_offset, model);
    
    const auto offsetted_truth_begin_itr = next(cbegin(truth), target_offset);
    
    const auto m1 = std::mismatch(cbegin(target), cend(target), offsetted_truth_begin_itr);
    
    if (m1.first == cend(target)) {
        return 0; // sequences are equal, can't do better than this
    }
    
    const auto m2 = std::mismatch(next(m1.first), cend(target), next(m1.second));
    
    if (m2.first == cend(target)) {
        // then there is only a single base difference between the sequences, can optimise
        const auto truth_index = distance(offsetted_truth_begin_itr, m1.second) + target_offset;
        
        if (truth_index < model.lhs_flank_size || truth_index >= (truth.size() - model.rhs_flank_size)) {
            return 0;
        }
        
        const auto target_index = distance(cbegin(target), m1.first);
        
        auto mispatch_penalty = target_qualities[target_index];
        
        if (model.snv_mask[truth_index] == *m1.first) {
            mispatch_penalty = std::min(target_qualities[target_index],
                                        static_cast<std::uint8_t>(model.snv_priors[truth_index]));
        }
        
        if (mispatch_penalty <= model.gap_open_penalties[truth_index]
            || !std::equal(next(m1.first), cend(target), m1.second)) {
            return Ln_probability[mispatch_penalty];
        }
        
        return Ln_probability[model.gap_open_penalties[truth_index]];
    }
    
    // TODO: we should be able to optimise the alignment based of the first mismatch postition
    
    return simd_align(truth, target, target_qualities, target_offset, model);
}

} // namespace hmm
} // namespace octopus
