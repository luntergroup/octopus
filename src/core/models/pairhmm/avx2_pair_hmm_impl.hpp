// Copyright (c) 2015-2019 Daniel Cooke and Gerton Lunter
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef avx2_pair_hmm_impl_hpp
#define avx2_pair_hmm_impl_hpp

#if __GNUC__ >= 6
    #pragma GCC diagnostic ignored "-Wignored-attributes"
#endif

#include <cstdint>
#include <immintrin.h>

#include <boost/container/small_vector.hpp>

#include "simd_pair_hmm.hpp"

namespace octopus { namespace hmm { namespace simd {

#ifdef __AVX2__

class InstructionSetPolicyAVX2
{
protected:
    using VectorType  = __m256i;
    using ScoreType   = short;
    using SmallVector = boost::container::small_vector<VectorType, 10000>;
    
    constexpr static int band_size = 16;
    constexpr static int trace_bits = 2;
    constexpr static ScoreType infinity = 0x7800;
    constexpr static ScoreType n_score = 2 << trace_bits;
    
    VectorType vectorise(ScoreType x) const noexcept
    {
        return _mm256_set1_epi16(penalty);
    }
    VectorType vectorise_zero_set_last(ScoreType x) const noexcept
    {
        return _mm256_set_epi16(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,penalty);
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
            default: return _extract<7>(a);
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
            default: return _insert<7>(a, i);
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
    template <int imm>
    VectorType _slli_si(const VectorType& a) const noexcept
    {
        return _mm256_slli_si256(a, imm);
    }
    template <int imm>
    VectorType _srli_si(const VectorType& a) const noexcept
    {
        return _mm256_srli_si256(a, imm);
    }
    template <int imm>
    VectorType _slli_epi(const VectorType& a) const noexcept
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

#endif /* __AVX2__ */

} // namespace simd
} // namespace hmm
} // namespace octopus

#endif
