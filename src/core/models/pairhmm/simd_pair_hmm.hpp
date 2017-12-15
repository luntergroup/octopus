// Copyright (c) 2017 Daniel Cooke and Gerton Lunter
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef simd_pair_hmm_hpp
#define simd_pair_hmm_hpp

#include <cstdint>

namespace octopus { namespace hmm { namespace simd {

constexpr int min_flank_pad() noexcept { return 8; }

int align(const char* truth, const char* target, const std::int8_t* qualities,
          int truth_len, int target_len,
          short gap_open, short gap_extend, short nuc_prior) noexcept;

int align(const char* truth, const char* target, const std::int8_t* qualities,
          int truth_len, int target_len,
          const std::int8_t* gap_open, short gap_extend, short nuc_prior) noexcept;

int align(const char* truth, const char* target, const std::int8_t* qualities,
          int truth_len, int target_len,
          const char* snv_mask, const std::int8_t* snv_prior,
          const std::int8_t* gap_open, short gap_extend, short nuc_prior) noexcept;

int align(const char* truth, const char* target, const std::int8_t* qualities,
          int truth_len, int target_len,
          const std::int8_t* gap_open, short gap_extend, short nuc_prior,
          int& first_pos, char* aln1, char* aln2) noexcept;

int align(const char* truth, const char* target, const std::int8_t* qualities,
          int truth_len, int target_len,
          const char* snv_mask, const std::int8_t* snv_prior,
          const std::int8_t* gap_open, short gap_extend, short nuc_prior,
          char* aln1, char* aln2, int& first_pos) noexcept;

int calculate_flank_score(int truth_len, int lhs_flank_len, int rhs_flank_len,
                          const char* target, const std::int8_t* quals,
                          const char* snv_mask, const std::int8_t* snv_prior,
                          const std::int8_t* gap_open, short gap_extend, short nuc_prior,
                          int first_pos, const char* aln1, const char* aln2,
                          int& target_mask_size) noexcept;

int calculate_flank_score(int truth_len, int lhs_flank_len, int rhs_flank_len,
                          const std::int8_t* quals, const std::int8_t* gap_open,
                          short gap_extend, short nuc_prior,
                          int first_pos, const char* aln1, const char* aln2,
                          int& target_mask_size) noexcept;

} // namespace simd
} // namespace hmm
} // namespace octopus

#endif
