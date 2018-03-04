// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

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

#include "utils/maths.hpp"
#include "simd_pair_hmm.hpp"

namespace octopus { namespace hmm {

using octopus::maths::constants::ln10Div10;

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
        const_cast<double&>(static_cast<std::array<double, N> const&>(result)[i]) = -ln10Div10<> * i;
    }
    return result;
}

template <typename T>
constexpr auto make_phred_to_ln_prob_lookup() noexcept
{
    return make_phred_to_ln_prob_lookup<num_values<T>()>();
}

bool target_overlaps_truth_flank(const std::string& truth, const std::string& target, const std::size_t target_offset,
                                 const MutationModel& model) noexcept
{
    constexpr auto pad = simd::min_flank_pad();
    return target_offset < (model.lhs_flank_size + pad)
           || (target_offset + target.size() + pad) > (truth.size() - model.rhs_flank_size);
}

bool use_adjusted_alignment_score(const std::string& truth, const std::string& target, const std::size_t target_offset,
                                  const MutationModel& model) noexcept
{
    return target_overlaps_truth_flank(truth, target, target_offset, model);
}

namespace debug {

void print_alignment(const std::vector<char>& align1, const std::vector<char>& align2)
{
    const auto isnt_null = [] (const char c) { return c != '\0'; };
    std::copy_if(std::cbegin(align1), std::cend(align1), std::ostreambuf_iterator<char> {std::cout}, isnt_null);
    std::cout << '\n';
    std::copy_if(std::cbegin(align2), std::cend(align2), std::ostreambuf_iterator<char> {std::cout}, isnt_null);
    std::cout << '\n';
}

} // namespace debug

auto make_cigar(const std::vector<char>& align1, const std::vector<char>& align2)
{
    assert(!align1.empty() && !align2.empty());
    auto align1_itr = std::cbegin(align1);
    auto align2_itr = std::cbegin(align2);
    const auto last_align1_itr = std::find_if_not(std::crbegin(align1), std::crend(align1), [] (auto x) { return x == 0; }).base();
    const auto last_align2_itr = std::next(align2_itr, std::distance(align1_itr, last_align1_itr));
    CigarString result {};
    result.reserve(std::distance(align1_itr, last_align1_itr));
    
    while (align1_itr != last_align1_itr) {
        const auto p = std::mismatch(align1_itr, last_align1_itr, align2_itr);
        if (p.first != align1_itr) {
            result.emplace_back(std::distance(align1_itr, p.first), CigarOperation::Flag::sequenceMatch);
            if (p.first == last_align1_itr) break;
        }
        const static auto is_gap = [] (const char b) noexcept { return b == '-'; };
        if (*p.first == '-') {
            const auto align1_gap_itr = std::find_if_not(std::next(p.first), last_align1_itr, is_gap);
            const auto gap_length = std::distance(p.first, align1_gap_itr);
            result.emplace_back(gap_length, CigarOperation::Flag::insertion);
            align1_itr = align1_gap_itr;
            align2_itr = std::next(p.second, gap_length);
        } else if (*p.second == '-') {
            const auto align2_gap_itr = std::find_if_not(std::next(p.second), last_align2_itr, is_gap);
            const auto gap_length = std::distance(p.second, align2_gap_itr);
            result.emplace_back(gap_length, CigarOperation::Flag::deletion);
            align1_itr = std::next(p.first, gap_length);
            align2_itr = align2_gap_itr;
        } else {
            const static auto is_mismatch = [] (char lhs, char rhs) noexcept { return lhs != rhs && lhs != '-' && rhs != '-'; };
            const auto p2 = std::mismatch(std::next(p.first), last_align1_itr, std::next(p.second), is_mismatch);
            result.emplace_back(std::distance(p.first, p2.first), CigarOperation::Flag::substitution);
            std::tie(align1_itr, align2_itr) = p2;
        }
    }
    return result;
}

auto simd_align(const std::string& truth, const std::string& target,
                const std::vector<std::uint8_t>& target_qualities,
                const std::size_t target_offset,
                const MutationModel& model) noexcept
{
    constexpr auto pad = simd::min_flank_pad();
    const auto truth_size  = static_cast<int>(truth.size());
    const auto target_size = static_cast<int>(target.size());
    const auto truth_alignment_size = static_cast<int>(target_size + 2 * pad - 1);
    const auto alignment_offset = std::max(0, static_cast<int>(target_offset) - pad);
    if (alignment_offset + truth_alignment_size > truth_size) {
        return std::numeric_limits<double>::lowest();
    }
    const auto qualities = reinterpret_cast<const std::int8_t*>(target_qualities.data());
    if (!use_adjusted_alignment_score(truth, target, target_offset, model)) {
        const auto score = simd::align(truth.data() + alignment_offset,
                                       target.data(),
                                       qualities,
                                       truth_alignment_size,
                                       target_size,
                                       model.snv_mask.data() + alignment_offset,
                                       model.snv_priors.data() + alignment_offset,
                                       model.gap_open.data() + alignment_offset,
                                       model.gap_extend, model.nuc_prior);
        return -ln10Div10<> * static_cast<double>(score);
    } else {
        thread_local std::vector<char> align1 {}, align2 {};
        const auto max_alignment_size = 2 * (target.size() + pad);
        align1.assign(max_alignment_size + 1, 0);
        align2.assign(max_alignment_size + 1, 0);
        int first_pos;
        const auto score = simd::align(truth.data() + alignment_offset,
                                       target.data(),
                                       qualities,
                                       truth_alignment_size,
                                       target_size,
                                       model.snv_mask.data() + alignment_offset,
                                       model.snv_priors.data() + alignment_offset,
                                       model.gap_open.data() + alignment_offset,
                                       model.gap_extend, model.nuc_prior,
                                       align1.data(), align2.data(), first_pos);
        auto lhs_flank_size = static_cast<int>(model.lhs_flank_size);
        if (lhs_flank_size < alignment_offset) {
            lhs_flank_size = 0;
        } else {
            lhs_flank_size -= alignment_offset;
            if (lhs_flank_size < 0) lhs_flank_size = 0;
        }
        auto rhs_flank_size = static_cast<int>(model.rhs_flank_size);
        if (alignment_offset + truth_alignment_size < truth_size - rhs_flank_size) {
            rhs_flank_size = 0;
        } else {
            rhs_flank_size += alignment_offset + truth_alignment_size;
            rhs_flank_size -= truth_size;
            if (rhs_flank_size < 0) rhs_flank_size = 0;
        }
        assert(lhs_flank_size >= 0 && rhs_flank_size >= 0);
        assert(align1.back() == 0); // required by calculate_flank_score
        int target_mask_size;
        auto flank_score = simd::calculate_flank_score(truth_alignment_size,
                                                       lhs_flank_size, rhs_flank_size,
                                                       target.data(), qualities,
                                                       model.snv_mask.data() + alignment_offset,
                                                       model.snv_priors.data() + alignment_offset,
                                                       model.gap_open.data() + alignment_offset,
                                                       model.gap_extend, model.nuc_prior,
                                                       first_pos,
                                                       align1.data(), align2.data(),
                                                       target_mask_size);
        const auto num_explained_bases = target_size - target_mask_size;
        constexpr int min_explained_bases {2};
        if (num_explained_bases < min_explained_bases) flank_score = 0;
        //assert(flank_score <= score);
        if (flank_score <= score) {
            return -ln10Div10<> * static_cast<double>(score - flank_score);
        } else {
            // Overflow has occurred when calculating score;
            return -ln10Div10<> * (flank_score + score);
        }
    }
}

auto simd_align_with_cigar(const std::string& truth, const std::string& target,
                           const std::vector<std::uint8_t>& target_qualities,
                           const std::size_t target_offset,
                           const MutationModel& model) noexcept
{
    constexpr auto pad = simd::min_flank_pad();
    const auto truth_size  = static_cast<int>(truth.size());
    const auto target_size = static_cast<int>(target.size());
    const auto truth_alignment_size = static_cast<int>(target_size + 2 * pad - 1);
    const auto alignment_offset = std::max(0, static_cast<int>(target_offset) - pad);
    if (alignment_offset + truth_alignment_size > truth_size) {
        return std::make_pair(CigarString {}, std::numeric_limits<double>::lowest());
    }
    const auto qualities = reinterpret_cast<const std::int8_t*>(target_qualities.data());
    thread_local std::vector<char> align1 {}, align2 {};
    const auto max_alignment_size = 2 * (target.size() + pad);
    align1.assign(max_alignment_size + 1, 0);
    align2.assign(max_alignment_size + 1, 0);
    int first_pos;
    auto score = simd::align(truth.data() + alignment_offset,
                             target.data(),
                             qualities,
                             truth_alignment_size,
                             target_size,
                             model.snv_mask.data() + alignment_offset,
                             model.snv_priors.data() + alignment_offset,
                             model.gap_open.data() + alignment_offset,
                             model.gap_extend, model.nuc_prior,
                             align1.data(), align2.data(), first_pos);
    if (use_adjusted_alignment_score(truth, target, target_offset, model)) {
        auto lhs_flank_size = static_cast<int>(model.lhs_flank_size);
        if (lhs_flank_size < alignment_offset) {
            lhs_flank_size = 0;
        } else {
            lhs_flank_size -= alignment_offset;
            if (lhs_flank_size < 0) lhs_flank_size = 0;
        }
        auto rhs_flank_size = static_cast<int>(model.rhs_flank_size);
        if (alignment_offset + truth_alignment_size < truth_size - rhs_flank_size) {
            rhs_flank_size = 0;
        } else {
            rhs_flank_size += alignment_offset + truth_alignment_size;
            rhs_flank_size -= truth_size;
            if (rhs_flank_size < 0) rhs_flank_size = 0;
        }
        assert(lhs_flank_size >= 0 && rhs_flank_size >= 0);
        assert(align1.back() == 0); // required by calculate_flank_score
        int target_mask_size;
        auto flank_score = simd::calculate_flank_score(truth_alignment_size,
                                                       lhs_flank_size, rhs_flank_size,
                                                       target.data(), qualities,
                                                       model.snv_mask.data() + alignment_offset,
                                                       model.snv_priors.data() + alignment_offset,
                                                       model.gap_open.data() + alignment_offset,
                                                       model.gap_extend, model.nuc_prior,
                                                       first_pos,
                                                       align1.data(), align2.data(),
                                                       target_mask_size);
        const auto num_explained_bases = target_size - target_mask_size;
        constexpr int min_explained_bases {2};
        if (num_explained_bases < min_explained_bases) flank_score = 0;
        if (flank_score <= score) {
            score -= flank_score;
        } else {
            score += flank_score;
            // Overflow has occurred when calculating score;
        }
    }
    return std::make_pair(make_cigar(align1, align2), -ln10Div10<> * static_cast<double>(score));
}

unsigned min_flank_pad() noexcept
{
    return simd::min_flank_pad();
}

void validate(const std::string& truth, const std::string& target,
              const std::vector<std::uint8_t>& target_qualities,
              const std::size_t target_offset,
              const MutationModel& model)
{
    if (target.size() != target_qualities.size()) {
        throw std::invalid_argument {"PairHMM::align: target size not equal to target base_qualities length"};
    }
    if (truth.size() != model.snv_priors.size()) {
        throw std::invalid_argument {"PairHMM::align: truth size not equal to snv priors length"};
    }
    if (truth.size() != model.gap_open.size()) {
        throw std::invalid_argument {"PairHMM::align: truth size not equal to gap open penalties length"};
    }
    if (target_offset + target.size() > truth.size()) {
        throw std::invalid_argument {"PairHMM::align: target is not contained by truth"};
    }
}

double evaluate(const std::string& target, const std::string& truth,
                const std::vector<std::uint8_t>& target_qualities,
                const std::size_t target_offset,
                const MutationModel& model)
{
    using std::cbegin; using std::cend; using std::next; using std::distance;
    static constexpr auto lnProbability = make_phred_to_ln_prob_lookup<std::uint8_t>();
    validate(truth, target, target_qualities, target_offset, model);
    const auto offsetted_truth_begin_itr = next(cbegin(truth), target_offset);
    const auto m1 = std::mismatch(cbegin(target), cend(target), offsetted_truth_begin_itr);
    if (m1.first == cend(target)) {
        return 0; // sequences are equal, can't do better than this
    }
    const auto m2 = std::mismatch(next(m1.first), cend(target), next(m1.second));
    if (m2.first == cend(target)) {
        // then there is only a single base difference between the sequences, can optimise
        const auto truth_mismatch_idx = distance(offsetted_truth_begin_itr, m1.second) + target_offset;
        if (truth_mismatch_idx < model.lhs_flank_size || truth_mismatch_idx >= (truth.size() - model.rhs_flank_size)) {
            return 0;
        }
        const auto target_index = distance(cbegin(target), m1.first);
        auto mispatch_penalty = target_qualities[target_index];
        if (model.snv_mask[truth_mismatch_idx] == *m1.first) {
            mispatch_penalty = std::min(target_qualities[target_index],
                                        static_cast<std::uint8_t>(model.snv_priors[truth_mismatch_idx]));
        }
        if (mispatch_penalty <= model.gap_open[truth_mismatch_idx]
            || !std::equal(next(m1.first), cend(target), m1.second)) {
            return lnProbability[mispatch_penalty];
        }
        return lnProbability[model.gap_open[truth_mismatch_idx]];
    }
    // TODO: we should be able to optimise the alignment based of the first mismatch postition
    return simd_align(truth, target, target_qualities, target_offset, model);
}

std::pair<CigarString, double>
align(const std::string& target, const std::string& truth,
      const std::vector<std::uint8_t>& target_qualities,
      std::size_t target_offset,
      const MutationModel& model)
{
    validate(truth, target, target_qualities, target_offset, model);
    if (std::equal(std::cbegin(target), std::cend(target), std::next(std::cbegin(truth), target_offset))) {
        return {{CigarOperation {static_cast<CigarOperation::Size>(target.size()), CigarOperation::Flag::sequenceMatch}}, 0};
    } else {
        return simd_align_with_cigar(truth, target, target_qualities, target_offset, model);
    }
}

double evaluate(const std::string& target, const std::string& truth, const VariableGapOpenMutationModel& model) noexcept
{
    assert(truth.size() == model.gap_open.size());
    using std::cbegin; using std::cend; using std::next; using std::distance;
    static constexpr auto lnProbability = make_phred_to_ln_prob_lookup<std::uint8_t>();
    const auto truth_begin = next(cbegin(truth), min_flank_pad());
    const auto m1 = std::mismatch(cbegin(target), cend(target), truth_begin);
    if (m1.first == cend(target)) {
        return 0; // sequences are equal, can't do better than this
    }
    const auto m2 = std::mismatch(next(m1.first), cend(target), next(m1.second));
    if (m2.first == cend(target)) {
        // then there is only a single base difference between the sequences, can optimise
        const auto truth_mismatch_idx = static_cast<std::size_t>(distance(cbegin(truth), m1.second));
        if (model.mutation <= model.gap_open[truth_mismatch_idx] || !std::equal(next(m1.first), cend(target), m1.second)) {
            return lnProbability[model.gap_open[truth_mismatch_idx]];
        }
        return lnProbability[model.mutation];
    }
    const auto truth_alignment_size = static_cast<int>(target.size() + 2 * min_flank_pad() - 1);
    thread_local std::vector<std::int8_t> dummy_qualities;
    dummy_qualities.assign(target.size(), model.mutation);
    auto score = simd::align(truth.c_str(), target.c_str(),
                             dummy_qualities.data(),
                             truth_alignment_size,
                             static_cast<int>(target.size()),
                             model.gap_open.data(),
                             model.gap_extend, 2);
    return -ln10Div10<> * static_cast<double>(score);
}

double evaluate(const std::string& target, const std::string& truth, const FlatGapMutationModel& model) noexcept
{
    using std::cbegin; using std::cend; using std::next; using std::distance;
    static constexpr auto lnProbability = make_phred_to_ln_prob_lookup<std::uint8_t>();
    const auto truth_begin = next(cbegin(truth), min_flank_pad());
    const auto m1 = std::mismatch(cbegin(target), cend(target), truth_begin);
    if (m1.first == cend(target)) {
        return 0; // sequences are equal, can't do better than this
    }
    const auto m2 = std::mismatch(next(m1.first), cend(target), next(m1.second));
    if (m2.first == cend(target)) {
        // then there is only a single base difference between the sequences, can optimise
        if (model.mutation <= model.gap_open || !std::equal(next(m1.first), cend(target), m1.second)) {
            return lnProbability[model.gap_open];
        }
        return lnProbability[model.mutation];
    }
    const auto truth_alignment_size = static_cast<int>(target.size() + 2 * min_flank_pad() - 1);
    thread_local std::vector<std::int8_t> dummy_qualities;
    dummy_qualities.assign(target.size(), model.mutation);
    auto score = simd::align(truth.c_str(), target.c_str(),
                             dummy_qualities.data(),
                             truth_alignment_size,
                             static_cast<int>(target.size()),
                             model.gap_open,
                             model.gap_extend, 2);
    return -ln10Div10<> * static_cast<double>(score);
}

} // namespace hmm
} // namespace octopus
