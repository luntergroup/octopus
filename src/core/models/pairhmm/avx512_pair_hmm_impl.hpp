// Copyright (c) 2015-2019 Daniel Cooke and Gerton Lunter
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef avx512_pair_hmm_impl_hpp
#define avx512_pair_hmm_impl_hpp

#if __GNUC__ >= 6
    #pragma GCC diagnostic ignored "-Wignored-attributes"
#endif

#include <cstdint>
#include <immintrin.h>

#include "simd_pair_hmm.hpp"

namespace octopus { namespace hmm { namespace simd {

namespace detail {

template <int n>
auto _left_shift_words(const __m512i& a) noexcept
{
    // TODO
    return a;
};

template <int n>
auto _right_shift_words(const __m512i& a) noexcept
{
    // TODO
    return a;
}

} // namespace detail

class AVX512PairHMMInstructionSet
{
protected:
    using VectorType = __m512i;
    using ScoreType  = short;
    
    constexpr static int word_size  = sizeof(ScoreType);
    constexpr static int band_size_ = sizeof(VectorType) / word_size;
    
    VectorType vectorise(ScoreType x) const noexcept
    {
        return _mm512_set1_epi16(x);
    }
    VectorType vectorise(const char* sequence) const noexcept
    {
        return _mm512_set_epi16(sequence[31], sequence[30], sequence[29], sequence[28],
                                sequence[27], sequence[26], sequence[25], sequence[24],
                                sequence[23], sequence[22], sequence[21], sequence[20],
                                sequence[19], sequence[18], sequence[17], sequence[16],
                                sequence[15], sequence[14], sequence[13], sequence[12],
                                sequence[11], sequence[10], sequence[9], sequence[8],
                                sequence[7], sequence[6], sequence[5], sequence[4],
                                sequence[3], sequence[2], sequence[1], sequence[0]);
    }
    VectorType vectorise_zero_set_last(ScoreType x) const noexcept
    {
        return _mm512_set_epi16(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,x);
    }
    template <int index>
    auto _extract(const VectorType a) const noexcept
    {
        static_assert(index < 32, "index must be less than 32");
        return _mm512_cvtsi512_si32(_mm512_srli_epi16(a, index)) & 0xFFFF;
    }
    auto _extract(const VectorType a, const int index) const noexcept
    {
        switch (index) {
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
            case 16:  return _extract<16>(a);
            case 17:  return _extract<17>(a);
            case 18:  return _extract<18>(a);
            case 19:  return _extract<19>(a);
            case 20:  return _extract<20>(a);
            case 21:  return _extract<21>(a);
            case 22:  return _extract<22>(a);
            case 23:  return _extract<23>(a);
            case 24:  return _extract<24>(a);
            case 25:  return _extract<25>(a);
            case 26:  return _extract<26>(a);
            case 27:  return _extract<27>(a);
            case 28:  return _extract<28>(a);
            case 29:  return _extract<29>(a);
            case 30:  return _extract<30>(a);
            case 31:  return _extract<31>(a);
            default: return _extract<31>(a);
        }
    }
    VectorType _insert_bottom(const VectorType& a, const ScoreType i) const noexcept
    {
        return _mm512_insert_epi16(a, i, 0);
    }
    VectorType _insert_top(const VectorType& a, const ScoreType i) const noexcept
    {
        return _mm512_insert_epi16(a, i, band_size_ - 1);
    }
    VectorType _add(const VectorType& lhs, const VectorType& rhs) const noexcept
    {
        return _mm512_add_epi16(lhs, rhs);
    }
    VectorType _and(const VectorType& lhs, const VectorType& rhs) const noexcept
    {
        return _mm512_and_si512(lhs, rhs);
    }
    VectorType _andnot(const VectorType& lhs, const VectorType& rhs) const noexcept
    {
        return _mm512_andnot_si512(lhs, rhs);
    }
    VectorType _or(const VectorType& lhs, const VectorType& rhs) const noexcept
    {
        return _mm512_or_si512(lhs, rhs);
    }
    VectorType _cmpeq(const VectorType& lhs, const VectorType& rhs) const noexcept
    {
        return _mm512_maskz_set1_epi16(_mm512_cmpeq_epi16_mask(lhs, rhs), 0xFFFF);
    }
    VectorType _left_shift_word(const VectorType& a) const noexcept
    {
        return detail::_left_shift_words<word_size>(a);
    }
    VectorType _right_shift_word(const VectorType& a) const noexcept
    {
        return detail::_right_shift_words<word_size>(a);
    }
    template <int n>
    VectorType _left_shift_bits(const VectorType& a) const noexcept
    {
        return _mm512_slli_epi16(a, n);
    }
    template <int n>
    VectorType _right_shift_bits(const VectorType& a) const noexcept
    {
        return _mm512_slli_epi16(a, n);
    }
    VectorType _min(const VectorType& lhs, const VectorType& rhs) const noexcept
    {
        return _mm512_min_epi16(lhs, rhs);
    }
    VectorType _max(const VectorType& lhs, const VectorType& rhs) const noexcept
    {
        return _mm512_max_epi16(lhs, rhs);
    }
};

using AVX512PairHMM = PairHMM<AVX512PairHMMInstructionSet>;

} // namespace simd
} // namespace hmm
} // namespace octopus

#endif
