/*****************************************************************************************************************
 This code is copyright (c) Gerton Lunter, Dec 2009
 It may not be distributed, made public, or used in other software without the permission of the copyright holder
 ******************************************************************************************************************/

#ifndef simd_pair_hmm_hpp
#define simd_pair_hmm_hpp

#include <cstdint>

namespace SimdPairHmm
{
    int align(const char* truth, const char* target, const std::int8_t* qualities,
              int truth_len, int target_len,
              const std::int8_t* gap_open, short gap_extend, short nuc_prior);
    
    int align(const char* truth, const char* target, const std::int8_t* qualities,
              int truth_len, int target_len,
              const std::int8_t* snv_prior, const std::int8_t* gap_open,
              short gap_extend, short nuc_prior);
    
    int align(const char* truth, const char* target, const std::int8_t* qualities,
              int truth_len, int target_len,
              const std::int8_t* gap_open, short gap_extend, short nuc_prior,
              char* aln1, char* aln2, int* first_pos);
    
    int align(const char* truth, const char* target, const std::int8_t* qualities,
              int truth_len, int target_len,
              const std::int8_t* snv_prior, const std::int8_t* gap_open,
              short gap_extend, short nuc_prior,
              char* aln1, char* aln2, int* first_pos);
    
    int calculate_flank_score(int truth_len, int lhs_flank_len, int rhs_flank_len,
                              const std::int8_t* quals, const std::int8_t* gap_open,
                              short gap_extend, short nuc_prior,
                              int first_pos, const char* aln1, const char* aln2);
} // namespace SimdPairHmm

#endif
