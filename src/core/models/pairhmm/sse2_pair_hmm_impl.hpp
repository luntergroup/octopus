// Copyright (c) 2015-2019 Daniel Cooke and Gerton Lunter
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef sse2_pair_hmm_impl_hpp
#define sse2_pair_hmm_impl_hpp

#if __GNUC__ >= 6
    #pragma GCC diagnostic ignored "-Wignored-attributes"
#endif

#include <cstdint>
#include <array>
#include <tuple>
#include <emmintrin.h>

#include "utils/array_tricks.hpp"
#include "simd_pair_hmm.hpp"

namespace octopus { namespace hmm { namespace simd {

template <unsigned MinBandSize>
class SSE2PairHMMInstructionSet
{
protected:
    using ScoreType = short;
    
    constexpr static int word_size = sizeof(ScoreType);

private:
    using BlockType = __m128i;
    
    constexpr static auto block_bytes_      = sizeof(BlockType);
    constexpr static auto block_words_      = block_bytes_ / word_size;
    constexpr static std::size_t num_blocks = MinBandSize / block_words_;
    
protected:
    using VectorType = std::array<BlockType, num_blocks>;
    
    constexpr static int band_size_ = sizeof(VectorType) / word_size;

private:
    constexpr static auto num_blocks_ = std::tuple_size<VectorType>::value;
    
    template <typename T, std::size_t...Is>
    VectorType vectorise_helper(const T* values, std::index_sequence<Is...>) const noexcept
    {
        return {{(static_cast<void>(Is),
            _mm_set_epi16(values[block_words_ * Is + 7],
                          values[block_words_ * Is + 6],
                          values[block_words_ * Is + 5],
                          values[block_words_ * Is + 4],
                          values[block_words_ * Is + 3],
                          values[block_words_ * Is + 2],
                          values[block_words_ * Is + 1],
                          values[block_words_ * Is + 0]
                          ))...}};
    }
    
protected:
    VectorType vectorise(ScoreType x) const noexcept
    {
        return make_array<num_blocks>(_mm_set1_epi16(x));
    }
    template <typename T>
    VectorType vectorise(const T* values) const noexcept
    {
        return vectorise_helper(values, std::make_index_sequence<num_blocks>());
    }
    VectorType vectorise_zero_set_last(ScoreType x) const noexcept
    {
        return make_array<num_blocks>(_mm_set_epi16(0,0,0,0,0,0,0,x), _mm_set_epi16(0,0,0,0,0,0,0,0));
    }
    template <int index>
    auto _extract(const VectorType a) const noexcept
    {
        constexpr static auto block_index = index / block_words_;
        static_assert(block_index < num_blocks_, "block index out range");
        constexpr static auto block_byte_index = index % block_words_;
        static_assert(block_byte_index < block_words_, "block byte index out range");
        return _mm_extract_epi16(std::get<block_index>(a), block_byte_index);
    }
    template <std::size_t... Is>
    auto _extract_helper(const VectorType& a, const int index, std::index_sequence<Is...>) const noexcept
    {
        if (index < band_size_ && index >= 0) {
            decltype(_extract<0>(a)) r;
            int unused[] = {(index == Is ? (r = _extract<Is>(a), 0) : 0)...};
            (void) unused;
            return r;
        } else {
            return _extract<band_size_ - 1>(a);
        }
    }
    auto _extract(const VectorType& a, const int index) const noexcept
    {
        return _extract_helper(a, index, std::make_index_sequence<band_size_> {});
    }
    VectorType _insert_bottom(VectorType a, const ScoreType i) const noexcept
    {
        std::get<0>(a) = _mm_insert_epi16(std::get<0>(a), i, 0);
        return a;
    }
    VectorType _insert_top(VectorType a, const ScoreType i) const noexcept
    {
        std::get<num_blocks_ - 1>(a) = _mm_insert_epi16(std::get<num_blocks_ - 1>(a), i, block_words_ - 1);
        return a;
    }
    VectorType _add(const VectorType& lhs, const VectorType& rhs) const noexcept
    {
        return transform([] (const auto& lhs, const auto& rhs) noexcept { return _mm_add_epi16(lhs, rhs); }, lhs, rhs);
    }
    VectorType _and(const VectorType& lhs, const VectorType& rhs) const noexcept
    {
        return transform([] (const auto& lhs, const auto& rhs) noexcept { return _mm_and_si128(lhs, rhs); }, lhs, rhs);
    }
    VectorType _andnot(const VectorType& lhs, const VectorType& rhs) const noexcept
    {
        return transform([] (const auto& lhs, const auto& rhs) noexcept { return _mm_andnot_si128(lhs, rhs); }, lhs, rhs);
    }
    VectorType _or(const VectorType& lhs, const VectorType& rhs) const noexcept
    {
        return transform([] (const auto& lhs, const auto& rhs) noexcept { return _mm_or_si128(lhs, rhs); }, lhs, rhs);
    }
    VectorType _cmpeq(const VectorType& lhs, const VectorType& rhs) const noexcept
    {
        return transform([] (const auto& lhs, const auto& rhs) noexcept { return _mm_cmpeq_epi16(lhs, rhs); }, lhs, rhs);
    }
    VectorType _left_shift_word(VectorType a) const noexcept
    {
        adjacent_apply_reverse([] (const auto& lhs, const auto& rhs) noexcept {
            return _mm_or_si128(_mm_slli_si128(rhs, word_size), _mm_srli_si128(lhs, block_bytes_ - word_size));
        }, a, a);
        std::get<0>(a) = _mm_slli_si128(std::get<0>(a), word_size);
        return a;
    }
    VectorType _right_shift_word(VectorType a) const noexcept
    {
        adjacent_apply([] (const auto& lhs, const auto& rhs) noexcept {
            return _mm_or_si128(_mm_srli_si128(lhs, word_size), _mm_slli_si128(rhs, block_bytes_ - word_size));
        }, a, a);
        std::get<num_blocks - 1>(a) = _mm_srli_si128(std::get<num_blocks - 1>(a), word_size);
        return a;
    }
    template <int n>
    VectorType _left_shift_bits(const VectorType& a) const noexcept
    {
        return transform([] (const auto& x) noexcept { return _mm_slli_epi16(x, n); }, a);
    }
    template <int n>
    VectorType _right_shift_bits(const VectorType& a) const noexcept
    {
        return transform([] (const auto& x) noexcept { return _mm_srli_epi16(x, n); }, a);
    }
    VectorType _min(const VectorType& lhs, const VectorType& rhs) const noexcept
    {
        return transform([] (const auto& lhs, const auto& rhs) noexcept { return _mm_min_epi16(lhs, rhs); }, lhs, rhs);
    }
    VectorType _max(const VectorType& lhs, const VectorType& rhs) const noexcept
    {
        return transform([] (const auto& lhs, const auto& rhs) noexcept { return _mm_max_epi16(lhs, rhs); }, lhs, rhs);
    }
};

template <unsigned MinBandSize>
using SSE2PairHMM = PairHMM<SSE2PairHMMInstructionSet<MinBandSize>>;

using FastestSSE2PairHMM = SSE2PairHMM<8>;

} // namespace simd
} // namespace hmm
} // namespace octopus

#endif
