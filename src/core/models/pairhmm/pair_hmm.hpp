// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef pair_hmm_hpp
#define pair_hmm_hpp

#include <string>
#include <vector>
#include <cstddef>
#include <cstdint>
#include <utility>
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

#include "basics/cigar_string.hpp"
#include "exceptions/program_error.hpp"
#include "utils/maths.hpp"
#include "simd_pair_hmm_factory.hpp"
#include "simd_pair_hmm_wrapper.hpp"

namespace octopus { namespace hmm {

// This is the minimum number of bases the truth must exceed the target either side around
// the mapped position
template <typename PairHMM>
unsigned min_flank_pad(const PairHMM& hmm) noexcept
{
    return hmm.band_size();
}

struct Alignment
{
    std::size_t target_offset;
    CigarString cigar;
    double likelihood;
};

class HMMOverflow : public ProgramError
{
public:
    using Sequence = std::string;
    
    HMMOverflow() = delete;
    HMMOverflow(const Sequence& target, const Sequence& truth) : target_ {target}, truth_ {truth} {}
    virtual ~HMMOverflow() override = default;
    
    const Sequence& target() const noexcept;
    const Sequence& truth() const noexcept;

private:
    const Sequence& target_, truth_;
    
    std::string do_why() const override { return "Pair HMM alignment overflowed"; }
    std::string do_where() const override { return "hmm::align"; }
};

struct NullType {};

using Penalty          = std::int8_t;
using PenaltyVector    = std::vector<Penalty>;
using NucleotideVector = std::vector<char>;

template <typename PenaltyOrPenaltyVector1,
          typename PenaltyOrPenaltyVector2,
          typename NucleotideVectorOrNull = NullType,
          typename PenaltyVectorOrNull = NullType,
          typename PenaltyOrNull = NullType,
          typename SizetOrNull = NullType>
struct Parameters
{
    PenaltyOrPenaltyVector1 gap_open;
    PenaltyOrPenaltyVector2 gap_extend;
    NucleotideVectorOrNull snv_mask = {};
    PenaltyVectorOrNull snv_priors = {};
    PenaltyOrNull mismatch = {};
    SizetOrNull lhs_flank_size = {}, rhs_flank_size = {};
    short nuc_prior = 2;
};

using MutationModel                  = Parameters<const PenaltyVector&, const PenaltyVector&, const NucleotideVector&, const PenaltyVector&, NullType, std::size_t>;
using VariableGapExtendMutationModel = Parameters<const PenaltyVector&, const PenaltyVector&, NullType, NullType, Penalty>;
using VariableGapOpenMutationModel   = Parameters<const PenaltyVector&, Penalty, NullType, NullType, Penalty>;
using FlatGapMutationModel           = Parameters<Penalty, Penalty, NullType, NullType, Penalty>;

using octopus::maths::constants::ln10Div10;

const static auto default_hmm = simd::make_simd_pair_hmm<16>();

namespace detail {

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

template <typename HMM, typename PairHMMParameters>
bool target_overlaps_truth_flank(const std::string& truth, const std::string& target, const std::size_t target_offset,
                                 const HMM& hmm, const PairHMMParameters& hmm_params) noexcept
{
    const auto pad = hmm.band_size();
    return target_offset < (hmm_params.lhs_flank_size + pad)
           || (target_offset + target.size() + pad) > (truth.size() - hmm_params.rhs_flank_size);
}

template <typename HMM, typename PairHMMParameters>
bool use_adjusted_alignment_score(const std::string& truth, const std::string& target, const std::size_t target_offset,
                                  const HMM& hmm, const PairHMMParameters& hmm_params) noexcept
{
    return target_overlaps_truth_flank(truth, target, target_offset, hmm, hmm_params);
}

namespace debug {

inline void print_alignment(const std::vector<char>& align1, const std::vector<char>& align2)
{
    const auto isnt_null = [] (const char c) { return c != '\0'; };
    std::copy_if(std::cbegin(align1), std::cend(align1), std::ostreambuf_iterator<char> {std::cout}, isnt_null);
    std::cout << '\n';
    std::copy_if(std::cbegin(align2), std::cend(align2), std::ostreambuf_iterator<char> {std::cout}, isnt_null);
    std::cout << '\n';
}

} // namespace debug

inline CigarString& make_cigar(const std::vector<char>& align1, const std::vector<char>& align2, CigarString& result)
{
    assert(!align1.empty() && !align2.empty());
    auto align1_itr = std::cbegin(align1);
    auto align2_itr = std::cbegin(align2);
    const auto last_align1_itr = std::find_if_not(std::crbegin(align1), std::crend(align1), [] (auto x) { return x == 0; }).base();
    const auto last_align2_itr = std::next(align2_itr, std::distance(align1_itr, last_align1_itr));
    result.clear();
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

inline CigarString make_cigar(const std::vector<char>& align1, const std::vector<char>& align2)
{
    CigarString result {};
    make_cigar(align1, align2, result);
    return result;
}

template <typename PairHMMParameters>
bool
is_in_flank(const std::size_t truth_idx,
            const std::size_t truth_length,
            const PairHMMParameters& model,
            std::true_type) noexcept
{
    return false;
}
template <typename PairHMMParameters>
bool
is_in_flank(const std::size_t truth_idx,
            const std::size_t truth_length,
            const PairHMMParameters& model,
            std::false_type) noexcept
{
    return truth_idx < model.lhs_flank_size || truth_idx >= (truth_length - model.rhs_flank_size);
}
template <typename PairHMMParameters>
bool
is_in_flank(const std::size_t truth_idx,
            const std::size_t truth_length,
            const PairHMMParameters& model) noexcept
{
    return is_in_flank(truth_idx, truth_length, model, std::is_same<decltype(model.lhs_flank_size), NullType> {});
}

template <typename Range>
auto get(const Range& values, std::size_t index, std::true_type) noexcept
{
    return values[index];
}
template <typename T>
auto get(const T& value, std::size_t, std::false_type) noexcept
{
    return value;
}
template <typename RangeOrConstant>
auto get(const RangeOrConstant& value, std::size_t index) noexcept
{
    return get(value, index, std::is_class<RangeOrConstant> {});
}

template <typename PairHMMParameters>
auto
get_mismatch_penalty(const char mismatched_target_base,
                     const std::uint8_t mismatched_target_base_quality,
                     std::size_t truth_mismatch_idx,
                     const PairHMMParameters& model,
                     std::true_type) noexcept
{
    return mismatched_target_base_quality;
}
template <typename PairHMMParameters>
auto
get_mismatch_penalty(const char mismatched_target_base,
                     const std::uint8_t mismatched_target_base_quality,
                     std::size_t truth_mismatch_idx,
                     const PairHMMParameters& model,
                     std::false_type) noexcept
{
    if (get(model.snv_mask, truth_mismatch_idx) == mismatched_target_base) {
        return std::min(mismatched_target_base_quality, static_cast<std::uint8_t>(get(model.snv_priors, truth_mismatch_idx)));
    } else {
        return mismatched_target_base_quality;
    }
}
template <typename PairHMMParameters>
auto
get_mismatch_penalty(const char mismatched_target_base,
                     const std::uint8_t mismatched_target_base_quality,
                     std::size_t truth_mismatch_idx,
                     const PairHMMParameters& model) noexcept
{
    return get_mismatch_penalty(mismatched_target_base, mismatched_target_base_quality, truth_mismatch_idx, model,
                                std::is_same<decltype(model.snv_mask), NullType> {});
}

template <typename Sequence1,
          typename Sequence2,
          typename PairHMMParameters>
std::pair<double, bool>
try_naive_evaluate(const Sequence1& truth,
                   const Sequence2& target,
                   const std::vector<std::uint8_t>& target_base_qualities,
                   const std::size_t target_offset,
                   const PairHMMParameters& model) noexcept
{
    using std::cbegin; using std::cend; using std::next; using std::distance;
    static constexpr auto ln_probability_table = make_phred_to_ln_prob_lookup<std::uint8_t>();
    const auto offsetted_truth_begin_itr = next(cbegin(truth), target_offset);
    const auto m1 = std::mismatch(cbegin(target), cend(target), offsetted_truth_begin_itr);
    if (m1.first == cend(target)) {
        return {0, true}; // sequences are equal, can't do better than this
    }
    const auto m2 = std::mismatch(next(m1.first), cend(target), next(m1.second));
    if (m2.first == cend(target)) {
        // then there is only a single base difference between the sequences, can optimise
        // target: ACGTACGT
        // truth:  ACGTTCGT
        const auto truth_mismatch_idx = distance(offsetted_truth_begin_itr, m1.second) + target_offset;
        if (is_in_flank(truth_mismatch_idx, truth.size(), model)) return {0, true};
        const auto target_mismatch_idx = distance(cbegin(target), m1.first);
        const auto mismatch_penalty = get_mismatch_penalty(*m1.first, target_base_qualities[target_mismatch_idx], truth_mismatch_idx, model);
        const auto gap_open_penalty = get(model.gap_open, truth_mismatch_idx);
        if (mismatch_penalty <= gap_open_penalty) {
            return {ln_probability_table[mismatch_penalty], true};
        } else {
            if (std::equal(next(m1.first), cend(target), m1.second)) {
                // target: AAAAGGGG
                // truth:  AAA GGGGG
                return {ln_probability_table[gap_open_penalty], true};
            } else if (std::equal(m1.first, cend(target), next(m1.second))) {
                // target: AAA GGGGG
                // truth:  AAAAGGGGG
                return {ln_probability_table[gap_open_penalty], true};
            } else if (mismatch_penalty <= (gap_open_penalty + get(model.gap_extend, truth_mismatch_idx))) {
                return {ln_probability_table[mismatch_penalty], true};
            }
        }
    }
    return {0, false};
}

template <typename Sequence1,
          typename Sequence2,
          typename PairHMMParameters>
bool
try_naive_align(const Sequence1& truth,
                const Sequence2& target,
                const std::vector<std::uint8_t>& target_base_qualities,
                const std::size_t target_offset,
                const PairHMMParameters& model,
                Alignment& result) noexcept
{
    if (std::equal(std::cbegin(target), std::cend(target), std::next(std::cbegin(truth), target_offset))) {
        result.likelihood = 0;
        result.target_offset = target_offset;
        result.cigar.assign({CigarOperation {static_cast<CigarOperation::Size>(target.size()), CigarOperation::Flag::sequenceMatch}});
        return true;
    } else {
        return false;
    }
}

template <typename Range>
auto data(const Range& values, std::size_t index, std::true_type) noexcept
{
    return values.data() + index;
}
template <typename T>
auto data(const T& value, std::size_t, std::false_type) noexcept
{
    return value;
}
template <typename RangeOrConstant>
auto data(const RangeOrConstant& value, std::size_t index) noexcept
{
    return data(value, index, std::is_class<RangeOrConstant> {});
}

template <typename Sequence1,
          typename Sequence2,
          typename PairHMM,
          typename PairHMMParameters>
auto
align(const Sequence1& truth,
      const Sequence2& target,
      const std::vector<std::uint8_t>& target_base_qualities,
      const int alignment_offset,
      const PairHMM& hmm,
      const PairHMMParameters& hmm_params,
      std::true_type) noexcept
{
    const auto truth_alignment_size = static_cast<int>(target.size() + 2 * hmm.band_size() - 1);
    const auto qualities = reinterpret_cast<const std::int8_t*>(target_base_qualities.data());
    return hmm.align(truth.data() + alignment_offset,
                     target.data(),
                     qualities,
                     truth_alignment_size,
                     target.size(),
                     data(hmm_params.gap_open, alignment_offset),
                     data(hmm_params.gap_extend, alignment_offset),
                     hmm_params.nuc_prior);
}
template <typename Sequence1,
          typename Sequence2,
          typename PairHMM,
          typename PairHMMParameters>
auto
align(const Sequence1& truth,
      const Sequence2& target,
      const std::vector<std::uint8_t>& target_base_qualities,
      const int alignment_offset,
      const PairHMM& hmm,
      const PairHMMParameters& hmm_params,
      std::false_type) noexcept
{
    const auto truth_alignment_size = static_cast<int>(target.size() + 2 * hmm.band_size() - 1);
    const auto qualities = reinterpret_cast<const std::int8_t*>(target_base_qualities.data());
    return hmm.align(truth.data() + alignment_offset,
                     target.data(),
                     qualities,
                     truth_alignment_size,
                     target.size(),
                     data(hmm_params.snv_mask, alignment_offset),
                     data(hmm_params.snv_priors, alignment_offset),
                     data(hmm_params.gap_open, alignment_offset),
                     data(hmm_params.gap_extend, alignment_offset),
                     hmm_params.nuc_prior);
}
template <typename Sequence1,
          typename Sequence2,
          typename PairHMM,
          typename PairHMMParameters>
auto
align(const Sequence1& truth,
      const Sequence2& target,
      const std::vector<std::uint8_t>& target_base_qualities,
      const int alignment_offset,
      const PairHMM& hmm,
      const PairHMMParameters& hmm_params) noexcept
{
    return align(truth, target, target_base_qualities, alignment_offset, hmm, hmm_params,
                 std::is_same<decltype(hmm_params.snv_mask), NullType> {});
}
template <typename Sequence1,
          typename Sequence2,
          typename PairHMM,
          typename PairHMMParameters>
auto
align(const Sequence1& truth,
      const Sequence2& target,
      const std::vector<std::uint8_t>& target_base_qualities,
      const int alignment_offset,
      const PairHMM& hmm,
      const PairHMMParameters& hmm_params,
      int& first_pos,
      std::vector<char>& align1,
      std::vector<char>& align2,
      std::true_type) noexcept
{
    const auto truth_alignment_size = static_cast<int>(target.size() + 2 * hmm.band_size() - 1);
    const auto qualities = reinterpret_cast<const std::int8_t*>(target_base_qualities.data());
    return hmm.align(truth.data() + alignment_offset,
                     target.data(),
                     qualities,
                     truth_alignment_size,
                     target.size(),
                     data(hmm_params.gap_open, alignment_offset),
                     data(hmm_params.gap_extend, alignment_offset),
                     hmm_params.nuc_prior,
                     first_pos,
                     align1.data(),
                     align2.data());
}
template <typename Sequence1,
          typename Sequence2,
          typename PairHMM,
          typename PairHMMParameters>
auto
align(const Sequence1& truth,
      const Sequence2& target,
      const std::vector<std::uint8_t>& target_base_qualities,
      const int alignment_offset,
      const PairHMM& hmm,
      const PairHMMParameters& hmm_params,
      int& first_pos,
      std::vector<char>& align1,
      std::vector<char>& align2,
      std::false_type) noexcept
{
    const auto truth_alignment_size = static_cast<int>(target.size() + 2 * hmm.band_size() - 1);
    const auto qualities = reinterpret_cast<const std::int8_t*>(target_base_qualities.data());
    return hmm.align(truth.data() + alignment_offset,
                     target.data(),
                     qualities,
                     truth_alignment_size,
                     target.size(),
                     data(hmm_params.snv_mask, alignment_offset),
                     data(hmm_params.snv_priors, alignment_offset),
                     data(hmm_params.gap_open, alignment_offset),
                     data(hmm_params.gap_extend, alignment_offset),
                     hmm_params.nuc_prior,
                     first_pos,
                     align1.data(),
                     align2.data());
}
template <typename Sequence1,
          typename Sequence2,
          typename PairHMM,
          typename PairHMMParameters>
auto
align(const Sequence1& truth,
      const Sequence2& target,
      const std::vector<std::uint8_t>& target_base_qualities,
      const int alignment_offset,
      const PairHMM& hmm,
      const PairHMMParameters& hmm_params,
      int& first_pos,
      std::vector<char>& align1,
      std::vector<char>& align2) noexcept
{
    return align(truth, target, target_base_qualities, alignment_offset, hmm, hmm_params,
                 first_pos, align1, align2,
                 std::is_same<decltype(hmm_params.snv_mask), NullType> {});
}

template <typename Sequence1,
          typename Sequence2,
          typename PairHMM,
          typename PairHMMParameters>
auto
calculate_flank_score(const Sequence1& truth,
                      const Sequence2& target,
                      const std::vector<std::uint8_t>& target_base_qualities,
                      const int alignment_offset,
                      const PairHMM& hmm,
                      const PairHMMParameters& hmm_params,
                      const int first_pos,
                      const std::vector<char>& align1,
                      const std::vector<char>& align2,
                      int& target_mask_size,
                      std::true_type) noexcept
{
    const auto truth_alignment_size = static_cast<int>(target.size() + 2 * hmm.band_size() - 1);
    const auto qualities = reinterpret_cast<const std::int8_t*>(target_base_qualities.data());
    const auto truth_size  = static_cast<int>(truth.size());
    auto lhs_flank_size = static_cast<int>(hmm_params.lhs_flank_size);
    if (lhs_flank_size < alignment_offset) {
        lhs_flank_size = 0;
    } else {
        lhs_flank_size -= alignment_offset;
        if (lhs_flank_size < 0) lhs_flank_size = 0;
    }
    auto rhs_flank_size = static_cast<int>(hmm_params.rhs_flank_size);
    if (alignment_offset + truth_alignment_size < truth_size - rhs_flank_size) {
        rhs_flank_size = 0;
    } else {
        rhs_flank_size += alignment_offset + truth_alignment_size;
        rhs_flank_size -= truth_size;
        if (rhs_flank_size < 0) rhs_flank_size = 0;
    }
    assert(lhs_flank_size >= 0 && rhs_flank_size >= 0);
    assert(align1.back() == 0); // required by calculate_flank_score
    return hmm.calculate_flank_score(truth_alignment_size,
                                     lhs_flank_size,
                                     rhs_flank_size,
                                     target.data(), qualities,
                                     data(hmm_params.gap_open, alignment_offset),
                                     data(hmm_params.gap_extend, alignment_offset),
                                     hmm_params.nuc_prior,
                                     first_pos,
                                     align1.data(), align2.data(),
                                     target_mask_size);
}
template <typename Sequence1,
          typename Sequence2,
          typename PairHMM,
          typename PairHMMParameters>
auto
calculate_flank_score(const Sequence1& truth,
                      const Sequence2& target,
                      const std::vector<std::uint8_t>& target_base_qualities,
                      const int alignment_offset,
                      const PairHMM& hmm,
                      const PairHMMParameters& hmm_params,
                      const int first_pos,
                      const std::vector<char>& align1,
                      const std::vector<char>& align2,
                      int& target_mask_size,
                      std::false_type) noexcept
{
    const auto truth_alignment_size = static_cast<int>(target.size() + 2 * hmm.band_size() - 1);
    const auto qualities = reinterpret_cast<const std::int8_t*>(target_base_qualities.data());
    const auto truth_size  = static_cast<int>(truth.size());
    auto lhs_flank_size = static_cast<int>(hmm_params.lhs_flank_size);
    if (lhs_flank_size < alignment_offset) {
        lhs_flank_size = 0;
    } else {
        lhs_flank_size -= alignment_offset;
        if (lhs_flank_size < 0) lhs_flank_size = 0;
    }
    auto rhs_flank_size = static_cast<int>(hmm_params.rhs_flank_size);
    if (alignment_offset + truth_alignment_size < truth_size - rhs_flank_size) {
        rhs_flank_size = 0;
    } else {
        rhs_flank_size += alignment_offset + truth_alignment_size;
        rhs_flank_size -= truth_size;
        if (rhs_flank_size < 0) rhs_flank_size = 0;
    }
    assert(lhs_flank_size >= 0 && rhs_flank_size >= 0);
    assert(align1.back() == 0); // required by calculate_flank_score
    return hmm.calculate_flank_score(truth_alignment_size,
                                     lhs_flank_size,
                                     rhs_flank_size,
                                     target.data(), qualities,
                                     data(hmm_params.snv_mask, alignment_offset),
                                     data(hmm_params.snv_priors, alignment_offset),
                                     data(hmm_params.gap_open, alignment_offset),
                                     data(hmm_params.gap_extend, alignment_offset),
                                     hmm_params.nuc_prior,
                                     first_pos,
                                     align1.data(), align2.data(),
                                     target_mask_size);
}
template <typename Sequence1,
          typename Sequence2,
          typename PairHMM,
          typename PairHMMParameters>
auto
calculate_flank_score(const Sequence1& truth,
                      const Sequence2& target,
                      const std::vector<std::uint8_t>& target_base_qualities,
                      const int alignment_offset,
                      const PairHMM& hmm,
                      const PairHMMParameters& hmm_params,
                      const int first_pos,
                      const std::vector<char>& align1,
                      const std::vector<char>& align2,
                      int& target_mask_size) noexcept
{
    return calculate_flank_score(truth, target, target_base_qualities, alignment_offset, hmm, hmm_params,
                                 first_pos, align1, align2,
                                 target_mask_size,
                                 std::is_same<decltype(hmm_params.snv_mask), NullType> {});
}

template <typename Sequence1,
          typename Sequence2,
          typename PairHMM,
          typename PairHMMParameters>
void
discount_flank_score(int& score,
                     const Sequence1& truth,
                     const Sequence2& target,
                     const std::vector<std::uint8_t>& target_base_qualities,
                     const std::size_t target_offset,
                     const int first_pos,
                     const std::vector<char>& align1,
                     const std::vector<char>& align2,
                     const PairHMM& hmm,
                     const PairHMMParameters& hmm_params,
                     std::true_type)
{}
template <typename Sequence1,
          typename Sequence2,
          typename PairHMM,
          typename PairHMMParameters>
void
discount_flank_score(int& score,
                     const Sequence1& truth,
                     const Sequence2& target,
                     const std::vector<std::uint8_t>& target_base_qualities,
                     const std::size_t target_offset,
                     const int first_pos,
                     const std::vector<char>& align1,
                     const std::vector<char>& align2,
                     const PairHMM& hmm,
                     const PairHMMParameters& hmm_params,
                     std::false_type)
{
    if (use_adjusted_alignment_score(truth, target, target_offset, hmm, hmm_params)) {
        const auto alignment_offset = std::max(0, static_cast<int>(target_offset) - hmm.band_size());
        int target_mask_size;
        auto flank_score = calculate_flank_score(truth, target, target_base_qualities, alignment_offset, hmm, hmm_params,
                                                 first_pos, align1, align2, target_mask_size);
        const auto num_explained_bases = target.size() - target_mask_size;
        constexpr int min_explained_bases {2};
        if (num_explained_bases < min_explained_bases) flank_score = 0;
        if (flank_score <= score) {
            score -= flank_score;
        } else {
            score += flank_score; // Overflow has occurred when calculating score;
        }
    }
}
template <typename Sequence1,
          typename Sequence2,
          typename PairHMM,
          typename PairHMMParameters>
void
discount_flank_score(int& score,
                     const Sequence1& truth,
                     const Sequence2& target,
                     const std::vector<std::uint8_t>& target_base_qualities,
                     const std::size_t target_offset,
                     const int first_pos,
                     const std::vector<char>& align1,
                     const std::vector<char>& align2,
                     const PairHMM& hmm,
                     const PairHMMParameters& hmm_params)
{
    discount_flank_score(score, truth, target, target_base_qualities, target_offset, first_pos, align1, align2, hmm, hmm_params,
                         std::is_same<decltype(hmm_params.lhs_flank_size), NullType> {});
}

template <typename Sequence1,
          typename Sequence2,
          typename PairHMM,
          typename PairHMMParameters>
auto
simd_evaluate_helper(const Sequence1& truth,
                     const Sequence2& target,
                     const std::vector<std::uint8_t>& target_base_qualities,
                     const std::size_t target_offset,
                     const PairHMM& hmm,
                     const PairHMMParameters& hmm_params,
                     std::true_type) noexcept
{
    const auto pad = hmm.band_size();
    const auto truth_size  = static_cast<int>(truth.size());
    const auto target_size = static_cast<int>(target.size());
    const auto truth_alignment_size = static_cast<int>(target_size + 2 * pad - 1);
    const auto alignment_offset = std::max(0, static_cast<int>(target_offset) - pad);
    if (alignment_offset + truth_alignment_size > truth_size) {
        return std::numeric_limits<double>::lowest();
    }
    auto score = align(truth, target, target_base_qualities, alignment_offset, hmm, hmm_params);
    return -ln10Div10<> * static_cast<double>(score);
}
template <typename Sequence1,
          typename Sequence2,
          typename PairHMM,
          typename PairHMMParameters>
auto
simd_evaluate_helper(const Sequence1& truth,
                     const Sequence2& target,
                     const std::vector<std::uint8_t>& target_base_qualities,
                     const std::size_t target_offset,
                     const PairHMM& hmm,
                     const PairHMMParameters& hmm_params,
                     std::false_type) noexcept
{
    const auto pad = hmm.band_size();
    const auto truth_size  = static_cast<int>(truth.size());
    const auto target_size = static_cast<int>(target.size());
    const auto truth_alignment_size = static_cast<int>(target_size + 2 * pad - 1);
    const auto alignment_offset = std::max(0, static_cast<int>(target_offset) - pad);
    if (alignment_offset + truth_alignment_size > truth_size) {
        return std::numeric_limits<double>::lowest();
    }
    if (!use_adjusted_alignment_score(truth, target, target_offset, hmm, hmm_params)) {
        auto score = align(truth, target, target_base_qualities, alignment_offset, hmm, hmm_params);
        return -ln10Div10<> * static_cast<double>(score);
    } else {
        thread_local std::vector<char> align1 {}, align2 {};
        const auto max_alignment_size = 2 * (target.size() + pad);
        align1.assign(max_alignment_size + 1, 0);
        align2.assign(max_alignment_size + 1, 0);
        int first_pos;
        const auto score = align(truth, target, target_base_qualities, alignment_offset, hmm, hmm_params,
                                 first_pos, align1, align2);
        if (first_pos == -1) {
            return std::numeric_limits<double>::lowest(); // overflow
        }
        assert(align1.back() == 0); // required by calculate_flank_score
        int target_mask_size;
        auto flank_score = calculate_flank_score(truth, target, target_base_qualities, alignment_offset, hmm, hmm_params,
                                                 first_pos, align1, align2, target_mask_size);
        const auto num_explained_bases = target_size - target_mask_size;
        constexpr int min_explained_bases {2};
        if (num_explained_bases < min_explained_bases) flank_score = 0;
        if (flank_score <= score) {
            return -ln10Div10<> * static_cast<double>(score - flank_score);
        } else {
            return -ln10Div10<> * (flank_score + score); // Overflow has occurred when calculating score;
        }
    }
}

template <typename Sequence1,
          typename Sequence2,
          typename PairHMM,
          typename PairHMMParameters>
auto
simd_evaluate(const Sequence1& truth,
              const Sequence2& target,
              const std::vector<std::uint8_t>& target_base_qualities,
              const std::size_t target_offset,
              const PairHMM& hmm,
              const PairHMMParameters& hmm_params) noexcept
{
    return simd_evaluate_helper(truth, target,target_base_qualities, target_offset, hmm, hmm_params,
                                std::is_same<decltype(hmm_params.lhs_flank_size), NullType> {});
}

template <typename Sequence1,
          typename Sequence2,
          typename PairHMM,
          typename PairHMMParameters>
void
simd_align(const Sequence1& truth,
           const Sequence2& target,
           const std::vector<std::uint8_t>& target_base_qualities,
           const std::size_t target_offset,
           const PairHMM& hmm,
           const PairHMMParameters& hmm_params,
           Alignment& result)
{
    const auto pad = hmm.band_size();
    const auto truth_size  = static_cast<int>(truth.size());
    const auto target_size = static_cast<int>(target.size());
    const auto truth_alignment_size = static_cast<int>(target_size + 2 * pad - 1);
    const auto alignment_offset = std::max(0, static_cast<int>(target_offset) - pad);
    if (alignment_offset + truth_alignment_size > truth_size) {
        result.target_offset = 0;
        result.likelihood = std::numeric_limits<double>::lowest();
        result.cigar = {};
        return;
    }
    thread_local std::vector<char> align1 {}, align2 {};
    const auto max_alignment_size = 2 * (target.size() + pad);
    align1.assign(max_alignment_size + 1, 0);
    align2.assign(max_alignment_size + 1, 0);
    int first_pos;
    auto score = align(truth, target, target_base_qualities, alignment_offset, hmm, hmm_params,
                       first_pos, align1, align2);
    if (first_pos == -1) {
        throw HMMOverflow {target, truth};
    }
    discount_flank_score(score, truth, target, target_base_qualities, target_offset,
                         first_pos, align1, align2, hmm, hmm_params);
    result.target_offset = target_offset - pad + first_pos;
    result.likelihood = -ln10Div10<> * static_cast<double>(score);
    make_cigar(align1, align2, result.cigar);
}

} // namespace detail

template <typename Sequence1,
          typename Sequence2,
          typename PairHMM,
          typename PairHMMParameters>
double
evaluate(const Sequence1& truth,
         const Sequence2& target,
         const std::vector<std::uint8_t>& target_base_qualities,
         const std::size_t target_offset,
         const PairHMM& hmm,
         const PairHMMParameters& model_params) noexcept
{
    auto p = detail::try_naive_evaluate(truth, target, target_base_qualities, target_offset, model_params);
    return p.second ? p.first : detail::simd_evaluate(truth, target, target_base_qualities, target_offset, hmm, model_params);
}

template <typename Sequence1,
          typename Sequence2,
          typename PairHMM,
          typename PairHMMParameters>
double
evaluate(const Sequence1& truth,
         const Sequence2& target,
         const PairHMM& hmm,
         const PairHMMParameters& model_params) noexcept
{
    thread_local std::vector<std::uint8_t> target_base_qualities;
    target_base_qualities.assign(target.size(), model_params.mismatch);
    return evaluate(truth, target, target_base_qualities, hmm.band_size(), hmm, model_params);
}

template <typename Sequence1,
          typename Sequence2,
          typename PairHMM,
          typename PairHMMParameters>
void
align(const Sequence1& target,
      const Sequence2& truth,
      const std::vector<std::uint8_t>& target_base_qualities,
      const std::size_t target_offset,
      const PairHMM& hmm,
      const PairHMMParameters& model_params,
      Alignment& result) noexcept
{
    if (!detail::try_naive_align(truth, target, target_base_qualities, target_offset, model_params, result)) {
        detail::simd_align(truth, target, target_base_qualities, target_offset, hmm, model_params, result);
    }
}

template <typename Sequence1,
          typename Sequence2,
          typename PairHMM,
          typename PairHMMParameters>
void
align(const Sequence1& target,
      const Sequence2& truth,
      const PairHMM& hmm,
      const PairHMMParameters& model_params,
      Alignment& result) noexcept
{
    thread_local std::vector<std::uint8_t> target_base_qualities;
    target_base_qualities.assign(target.size(), model_params.mismatch);
    align(truth, target, target_base_qualities, hmm.band_size(), hmm, model_params, result);
}

template <typename Parameters,
          int BandSize = 0,
          typename ScoreType = NullType>
class PairHMM
{
public:
    using ParameterType = Parameters;
    
    PairHMM() = default;
    PairHMM(unsigned min_band_size) { reset(min_band_size); }
    PairHMM(const Parameters& params, unsigned min_band_size = 8)
    {
        reset(min_band_size);
        set(params);
    }
    
    PairHMM(const PairHMM&)            = default;
    PairHMM& operator=(const PairHMM&) = default;
    PairHMM(PairHMM&&)                 = default;
    PairHMM& operator=(PairHMM&&)      = default;
    
    ~PairHMM() = default;
    
    int band_size() const noexcept { return hmm_.band_size(); }
    
    void set(const Parameters& params) noexcept { params_ = std::addressof(params); }
    
    template <typename Sequence1,
              typename Sequence2>
    double
    evaluate(const Sequence1& target,
             const Sequence2& truth,
             const std::vector<std::uint8_t>& target_base_qualities,
             const std::size_t target_offset) const noexcept
    {
        assert(params_);
        return octopus::hmm::evaluate(target, truth, target_base_qualities, target_offset, hmm_, *params_);
    }
    
    template <typename Sequence1,
              typename Sequence2>
    double
    evaluate(const Sequence1& target,
             const Sequence2& truth) const noexcept
    {
        assert(params_);
        return octopus::hmm::evaluate(target, truth, hmm_, *params_);
    }
    
    template <typename Sequence1,
              typename Sequence2>
    void
    align(const Sequence1& target,
          const Sequence2& truth,
          const std::vector<std::uint8_t>& target_base_qualities,
          const std::size_t target_offset,
          Alignment& result) const
    {
        assert(params_);
        octopus::hmm::align(truth, target, target_base_qualities, target_offset, hmm_, *params_, result);
    }
    template <typename Sequence1,
              typename Sequence2>
    Alignment
    align(const Sequence1& target,
          const Sequence2& truth,
          const std::vector<std::uint8_t>& target_base_qualities,
          const std::size_t target_offset) const
    {
        Alignment result {};
        this->align(target, truth, target_base_qualities, target_offset, result);
        return result;
    }
    
    template <typename Sequence1,
              typename Sequence2>
    void
    align(const Sequence1& target,
          const Sequence2& truth,
          Alignment& result) const
    {
        assert(params_);
        octopus::hmm::align(truth, target, hmm_, *params_, result);
    }
    template <typename Sequence1,
              typename Sequence2>
    Alignment
    align(const Sequence1& target,
          const Sequence2& truth) const
    {
        Alignment result {};
        this->align(target, truth, result);
        return result;
    }
    
private:
    using SIMDHMM = std::conditional_t<(BandSize > 0), simd::SimdPairHMM<BandSize, ScoreType>, simd::PairHMMWrapper>;
    
    static constexpr bool is_static = BandSize > 0 && !std::is_same<ScoreType, NullType>::value;
    
    SIMDHMM hmm_;
    const Parameters* params_ = nullptr;
    
    void reset(unsigned min_band_size) { reset(min_band_size, std::conditional_t<is_static, std::true_type, std::false_type> {}); }
    void reset(unsigned min_band_size, std::true_type) const noexcept {}
    void reset(unsigned min_band_size, std::false_type) { hmm_.reset(min_band_size, score_precision()); }
    simd::PairHMMWrapper::ScorePrecision score_precision(NullType) const noexcept
    {
        return simd::PairHMMWrapper::ScorePrecision::int16;
    }
    simd::PairHMMWrapper::ScorePrecision score_precision(short) const noexcept
    {
        return simd::PairHMMWrapper::ScorePrecision::int16;
    }
    simd::PairHMMWrapper::ScorePrecision score_precision(int) const noexcept
    {
        return simd::PairHMMWrapper::ScorePrecision::int32;
    }
    simd::PairHMMWrapper::ScorePrecision score_precision() const noexcept
    {
        return score_precision(ScoreType {});
    }
};

} // namespace hmm
} // namespace octopus

#endif
