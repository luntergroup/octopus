// Copyright (c) 2015-2019 Daniel Cooke and Gerton Lunter
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef avx2_pair_hmm_impl_hpp
#define avx2_pair_hmm_impl_hpp

#if __GNUC__ >= 6
#pragma GCC diagnostic ignored "-Wignored-attributes"
#endif

#include <cstdint>
#include <immintrin.h>

#include "simd_pair_hmm.hpp"

namespace octopus { namespace hmm { namespace simd {

namespace detail {

template <int index>
auto _left_shift(const __m256i& a) noexcept
{
    static_assert(index < 16, "index must be less than 16");
    return _mm256_alignr_epi8(a, _mm256_permute2x128_si256(a, a, _MM_SHUFFLE(0, 0, 2, 0)), 16 - index);
};
template <>
auto _left_shift<16>(const __m256i& a) noexcept
{
    return _mm256_permute2x128_si256(a, a, _MM_SHUFFLE(0, 0, 2, 0));
}

template <int index>
auto _right_shift(const __m256i& a) noexcept
{
    static_assert(index < 16, "index must be less than 16");
    return _mm256_alignr_epi8(_mm256_permute2x128_si256(a, a, _MM_SHUFFLE(2, 0, 0, 1)), a, index);
}
template <>
auto _right_shift<16>(const __m256i& a) noexcept
{
    return _mm256_permute2x128_si256(a, a, _MM_SHUFFLE(2, 0, 0, 1));
}

} // namespace detail

class InstructionSetPolicyAVX2
{
protected:
    using VectorType = __m256i;
    using ScoreType  = short;
    
    VectorType vectorise(ScoreType x) const noexcept
    {
        return _mm256_set1_epi16(x);
    }
    VectorType vectorise_zero_set_last(ScoreType x) const noexcept
    {
        return _mm256_set_epi16(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,x);
    }
    VectorType vectorise_reverse(const char* sequence) const noexcept
    {
        return _mm256_set_epi16(sequence[15], sequence[14], sequence[13], sequence[12],
                                sequence[11], sequence[10], sequence[9], sequence[8],
                                sequence[7], sequence[6], sequence[5], sequence[4],
                                sequence[3], sequence[2], sequence[1], sequence[0]);
    }
    VectorType vectorise_reverse_lshift(const std::int8_t* values, const int shift) const noexcept
    {
        return _mm256_set_epi16(values[15] << shift, values[14] << shift, values[13] << shift, values[12] << shift,
                                values[11] << shift, values[10] << shift, values[9] << shift, values[8] << shift,
                                values[7] << shift, values[6] << shift, values[5] << shift, values[4] << shift,
                                values[3] << shift, values[2] << shift, values[1] << shift, values[0] << shift);
    }
    VectorType vectorise_reverse_lshift(const std::int8_t value, const int shift) const noexcept
    {
        return vectorise(value << shift);
    }
    template <int imm>
    auto _extract(const VectorType a) const noexcept
    {
        return _mm256_extract_epi16(a, imm);
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
            case 8:  return _extract<8>(a);
            case 9:  return _extract<9>(a);
            case 10:  return _extract<10>(a);
            case 11:  return _extract<11>(a);
            case 12:  return _extract<12>(a);
            case 13:  return _extract<13>(a);
            case 14:  return _extract<14>(a);
            case 15:  return _extract<15>(a);
            default: return _extract<15>(a);
        }
    }
    template <int imm, typename T>
    VectorType _insert(const VectorType& a, T i) const noexcept
    {
        return _mm256_insert_epi16(a, i, imm);
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
            case 8:  return _insert<8>(a, i);
            case 9:  return _insert<9>(a, i);
            case 10:  return _insert<10>(a, i);
            case 11:  return _insert<11>(a, i);
            case 12:  return _insert<12>(a, i);
            case 13:  return _insert<13>(a, i);
            case 14:  return _insert<14>(a, i);
            case 15:  return _insert<15>(a, i);
            default: return _insert<15>(a, i);
        }
    }
    VectorType _add(const VectorType& lhs, const VectorType& rhs) const noexcept
    {
        return _mm256_add_epi16(lhs, rhs);
    }
    VectorType _and(const VectorType& lhs, const VectorType& rhs) const noexcept
    {
        return _mm256_and_si256(lhs, rhs);
    }
    VectorType _andnot(const VectorType& lhs, const VectorType& rhs) const noexcept
    {
        return _mm256_andnot_si256(lhs, rhs);
    }
    VectorType _or(const VectorType& lhs, const VectorType& rhs) const noexcept
    {
        return _mm256_or_si256(lhs, rhs);
    }
    VectorType _cmpeq(const VectorType& lhs, const VectorType& rhs) const noexcept
    {
        return _mm256_cmpeq_epi16(lhs, rhs);
    }
    template <int index>
    VectorType _left_shift(const VectorType& a) const noexcept
    {
        return detail::_left_shift<index>(a);
    }
    template <int index>
    VectorType _right_shift(const VectorType& a) const noexcept
    {
        return detail::_right_shift<index>(a);
    }
    template <int imm>
    VectorType _left_shift_words(const VectorType& a) const noexcept
    {
        return _mm256_slli_epi16(a, imm);
    }
    VectorType _min(const VectorType& lhs, const VectorType& rhs) const noexcept
    {
        return _mm256_min_epi16(lhs, rhs);
    }
    VectorType _max(const VectorType& lhs, const VectorType& rhs) const noexcept
    {
        return _mm256_max_epi16(lhs, rhs);
    }
};

using AVX2PairHMM = PairHMM<InstructionSetPolicyAVX2>;

} // namespace simd
} // namespace hmm
} // namespace octopus

#endif
