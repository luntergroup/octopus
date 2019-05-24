// Copyright (c) 2015-2019 Daniel Cooke and Gerton Lunter
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef sse2_pair_hmm_impl_hpp
#define sse2_pair_hmm_impl_hpp

#if __GNUC__ >= 6
    #pragma GCC diagnostic ignored "-Wignored-attributes"
#endif

#include <cstdint>
#include <emmintrin.h>

#include "simd_pair_hmm.hpp"

namespace octopus { namespace hmm { namespace simd {

class InstructionSetPolicySSE2
{
protected:
    using VectorType  = __m128i;
    using ScoreType   = short;
    
    VectorType vectorise(ScoreType x) const noexcept
    {
        return _mm_set1_epi16(x);
    }
    VectorType vectorise_zero_set_last(ScoreType x) const noexcept
    {
        return _mm_set_epi16(0,0,0,0,0,0,0,x);
    }
    VectorType vectorise_reverse(const char* sequence) const noexcept
    {
        return _mm_set_epi16(sequence[7], sequence[6], sequence[5], sequence[4],
                             sequence[3], sequence[2], sequence[1], sequence[0]);
    }
    VectorType vectorise_reverse_lshift(const std::int8_t* values, const int shift) const noexcept
    {
        return _mm_set_epi16(values[7] << shift, values[6] << shift, values[5] << shift, values[4] << shift,
                             values[3] << shift, values[2] << shift, values[1] << shift, values[0] << shift);
    }
    VectorType vectorise_reverse_lshift(const std::int8_t value, const int shift) const noexcept
    {
        return vectorise(value << shift);
    }
    template <int imm>
    auto _extract(const VectorType a) const noexcept
    {
        return _mm_extract_epi16(a, imm);
    }
    auto _extract(const VectorType a, const int n) const noexcept
    {
        switch (n) {
            case 0:  return _extract<0>(a);
            case 1:  return _extract<1>(a);
            case 2:  return _extract<2>(a);
            case 3:  return _extract<3>(a);
            case 4:  return _extract<4>(a);
            case 5:  return _extract<5>(a);
            case 6:  return _extract<6>(a);
            case 7:  return _extract<7>(a);
            default: return _extract<7>(a);
        }
    }
    template <int imm, typename T>
    VectorType _insert(const VectorType& a, T i) const noexcept
    {
        return _mm_insert_epi16(a, i, imm);
    }
    template <typename T>
    VectorType _insert(const VectorType& a, const T i, const int n) const noexcept
    {
        switch (n) {
            case 0:  return _insert<0>(a, i);
            case 1:  return _insert<1>(a, i);
            case 2:  return _insert<2>(a, i);
            case 3:  return _insert<3>(a, i);
            case 4:  return _insert<4>(a, i);
            case 5:  return _insert<5>(a, i);
            case 6:  return _insert<6>(a, i);
            case 7:  return _insert<7>(a, i);
            default: return _insert<7>(a, i);
        }
    }
    VectorType _add(const VectorType& lhs, const VectorType& rhs) const noexcept
    {
        return _mm_add_epi16(lhs, rhs);
    }
    VectorType _and(const VectorType& lhs, const VectorType& rhs) const noexcept
    {
        return _mm_and_si128(lhs, rhs);
    }
    VectorType _andnot(const VectorType& lhs, const VectorType& rhs) const noexcept
    {
        return _mm_andnot_si128(lhs, rhs);
    }
    VectorType _or(const VectorType& lhs, const VectorType& rhs) const noexcept
    {
        return _mm_or_si128(lhs, rhs);
    }
    VectorType _cmpeq(const VectorType& lhs, const VectorType& rhs) const noexcept
    {
        return _mm_cmpeq_epi16(lhs, rhs);
    }
    template <int imm>
    VectorType _left_shift(const VectorType& a) const noexcept
    {
        return _mm_slli_si128(a, imm);
    }
    template <int imm>
    VectorType _right_shift(const VectorType& a) const noexcept
    {
        return _mm_srli_si128(a, imm);
    }
    template <int imm>
    VectorType _left_shift_words(const VectorType& a) const noexcept
    {
        return _mm_slli_epi16(a, imm);
    }
    VectorType _min(const VectorType& lhs, const VectorType& rhs) const noexcept
    {
        return _mm_min_epi16(lhs, rhs);
    }
    VectorType _max(const VectorType& lhs, const VectorType& rhs) const noexcept
    {
        return _mm_max_epi16(lhs, rhs);
    }
};

using SSE2PairHMM = PairHMM<InstructionSetPolicySSE2>;

} // namespace simd
} // namespace hmm
} // namespace octopus

#endif
