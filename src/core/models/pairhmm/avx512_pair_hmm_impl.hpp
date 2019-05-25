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

template <int index>
auto _left_shift(const __m512i& a) noexcept
{
    // TODO
    return a;
};

template <int index>
auto _right_shift(const __m512i& a) noexcept
{
    // TODO
    return a;
}

} // namespace detail

class InstructionSetPolicyAVX512
{
protected:
    using VectorType = __m512i;
    using ScoreType  = short;
    
    VectorType vectorise(ScoreType x) const noexcept
    {
        return _mm512_set1_epi16(x);
    }
    VectorType vectorise_zero_set_last(ScoreType x) const noexcept
    {
        return _mm512_set_epi16(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,x);
    }
    VectorType vectorise_reverse(const char* sequence) const noexcept
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
    VectorType vectorise_reverse_lshift(const std::int8_t* values, const int shift) const noexcept
    {
        return _mm512_set_epi16(values[31] << shift, values[30] << shift, values[29] << shift, values[28] << shift,
                                values[27] << shift, values[26] << shift, values[25] << shift, values[24] << shift,
                                values[23] << shift, values[22] << shift, values[21] << shift, values[20] << shift,
                                values[19] << shift, values[18] << shift, values[17] << shift, values[16] << shift,
                                values[15] << shift, values[14] << shift, values[13] << shift, values[12] << shift,
                                values[11] << shift, values[10] << shift, values[9] << shift, values[8] << shift,
                                values[7] << shift, values[6] << shift, values[5] << shift, values[4] << shift,
                                values[3] << shift, values[2] << shift, values[1] << shift, values[0] << shift);
    }
    VectorType vectorise_reverse_lshift(const std::int8_t value, const int shift) const noexcept
    {
        return vectorise(value << shift);
    }
    template <int index>
    auto _extract(const VectorType a) const noexcept
    {
        static_assert(index < 32, "index must be less than 32");
        // TODO
        //return _mm512_extract_epi16(a, index);
        return a;
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
    template <int index, typename T>
    VectorType _insert(const VectorType& a, T i) const noexcept
    {
        static_assert(index < 32, "index must be less than 32");
        return _mm512_insert_epi16(a, i, index);
    }
    template <typename T>
    VectorType _insert(const VectorType& a, const T i, const int index) const noexcept
    {
        switch (index) {
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
            case 16:  return _insert<16>(a, i);
            case 17:  return _insert<17>(a, i);
            case 18:  return _insert<18>(a, i);
            case 19:  return _insert<19>(a, i);
            case 20:  return _insert<20>(a, i);
            case 21:  return _insert<21>(a, i);
            case 22:  return _insert<22>(a, i);
            case 23:  return _insert<23>(a, i);
            case 24:  return _insert<24>(a, i);
            case 25:  return _insert<25>(a, i);
            case 26:  return _insert<26>(a, i);
            case 27:  return _insert<27>(a, i);
            case 28:  return _insert<28>(a, i);
            case 29:  return _insert<29>(a, i);
            case 30:  return _insert<30>(a, i);
            case 31:  return _insert<31>(a, i);
            default: return _insert<31>(a, i);
        }
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
        // TODO: use _mm512_cmpeq_epi16_mask?
        //return _mm512_cmpeq_epi16(lhs, rhs);
        return lhs;
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
    template <int index>
    VectorType _left_shift_words(const VectorType& a) const noexcept
    {
        return _mm512_slli_epi16(a, index);
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

using AVX512PairHMM = PairHMM<InstructionSetPolicyAVX512>;

} // namespace simd
} // namespace hmm
} // namespace octopus

#endif
