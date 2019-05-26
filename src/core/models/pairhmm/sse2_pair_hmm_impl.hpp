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
#include <type_traits>
#include <emmintrin.h>

#include "utils/array_tricks.hpp"

namespace octopus { namespace hmm { namespace simd {

template <typename T>
constexpr bool is_short_or_int = std::is_same<T, short>::value || std::is_same<T, int>::value;

template <unsigned MinBandSize,
          typename ScoreTp = short>
class SSE2PairHMMInstructionSet
{
protected:
    using ScoreType = ScoreTp;
    
    static_assert(is_short_or_int<ScoreType>, "ScoreType not short or int");
    
    constexpr static int word_size = sizeof(ScoreType);
    
    constexpr static const char* name = "SSE2";

private:
    using BlockType = __m128i;
    
    constexpr static auto block_bytes_      = sizeof(BlockType);
    constexpr static auto block_words_      = block_bytes_ / word_size;
    constexpr static std::size_t num_blocks = MinBandSize / block_words_;
    
protected:
    using VectorType = std::array<BlockType, num_blocks>;
    
    constexpr static int band_size = sizeof(VectorType) / word_size;

private:
    constexpr static auto num_blocks_ = std::tuple_size<VectorType>::value;
    
    static VectorType do_vectorise(ScoreType x, short) noexcept
    {
        return make_array<num_blocks>(_mm_set1_epi16(x));
    }
    static VectorType do_vectorise(ScoreType x, int) noexcept
    {
        return make_array<num_blocks>(_mm_set1_epi32(x));
    }
protected:
    static VectorType vectorise(ScoreType x) noexcept
    {
        return do_vectorise(x, ScoreType {});
    }
private:
    template <typename T, std::size_t...Is>
    static VectorType do_vectorise(const T* values, std::index_sequence<Is...>, short) noexcept
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
    template <typename T, std::size_t...Is>
    static VectorType do_vectorise(const T* values, std::index_sequence<Is...>, int) noexcept
    {
        return {{(static_cast<void>(Is),
        _mm_set_epi32(values[block_words_ * Is + 3],
                      values[block_words_ * Is + 2],
                      values[block_words_ * Is + 1],
                      values[block_words_ * Is + 0]
        ))...}};
    }
protected:
    template <typename T>
    static VectorType vectorise(const T* values) noexcept
    {
        return do_vectorise(values, std::make_index_sequence<num_blocks>(), ScoreType {});
    }
private:
    static VectorType do_vectorise_zero_set_last(ScoreType x, short) noexcept
    {
        return make_array<num_blocks>(_mm_set_epi16(0,0,0,0,0,0,0,x), _mm_set_epi16(0,0,0,0,0,0,0,0));
    }
    static VectorType do_vectorise_zero_set_last(ScoreType x, int) noexcept
    {
        return make_array<num_blocks>(_mm_set_epi32(0,0,0,x), _mm_set_epi32(0,0,0,0));
    }
protected:
    static VectorType vectorise_zero_set_last(ScoreType x) noexcept
    {
        return do_vectorise_zero_set_last(x, ScoreType {});
    }
private:
    template <int index>
    static auto do_extract(const BlockType& a, short) noexcept
    {
        return _mm_extract_epi16(a, index);
    }
    template <int index>
    static auto do_extract(const BlockType& a, int) noexcept
    {
        return _mm_extract_epi32(a, index);
    }
protected:
    template <int index>
    static auto _extract(const VectorType a) noexcept
    {
        constexpr static auto block_index = index / block_words_;
        static_assert(block_index < num_blocks_, "block index out range");
        constexpr static auto block_byte_index = index % block_words_;
        static_assert(block_byte_index < block_words_, "block byte index out range");
        return do_extract<block_byte_index>(std::get<block_index>(a), ScoreType {});
    }
private:
    template <std::size_t... Is>
    static auto _extract_helper(const VectorType& a, const int index, std::index_sequence<Is...>) noexcept
    {
        if (index < band_size && index >= 0) {
            decltype(_extract<0>(a)) r;
            int unused[] = {(index == Is ? (r = _extract<Is>(a), 0) : 0)...};
            (void) unused;
            return r;
        } else {
            return _extract<band_size - 1>(a);
        }
    }
protected:
    static auto _extract(const VectorType& a, const int index) noexcept
    {
        return _extract_helper(a, index, std::make_index_sequence<band_size> {});
    }
private:
    template <int index, typename T>
    static BlockType do_insert(const BlockType& a, T value, short) noexcept
    {
        return _mm_insert_epi16(a, value, index);
    }
    template <int index, typename T>
    static BlockType do_insert(const BlockType& a, T value, int) noexcept
    {
        return _mm_insert_epi32(a, value, index);
    }
protected:
    template <int index, typename T>
    static VectorType _insert(VectorType a, T value) noexcept
    {
        constexpr static auto block_index = index / block_words_;
        static_assert(block_index < num_blocks_, "block index out range");
        constexpr static auto block_byte_index = index % block_words_;
        static_assert(block_byte_index < block_words_, "block byte index out range");
        std::get<block_index>(a) = do_insert<block_byte_index>(std::get<block_index>(a), value, ScoreType {});
        return a;
    }
private:
    template <typename T, std::size_t... Is>
    static auto _insert_helper(const VectorType& a, const T& value, const int index, std::index_sequence<Is...>) noexcept
    {
        if (index < band_size && index >= 0) {
            decltype(_insert<0>(a, value)) r;
            int unused[] = {(index == Is ? (r = _insert<Is>(a, value), 0) : 0)...};
            (void) unused;
            return r;
        } else {
            return _insert<band_size - 1>(a, value);
        }
    }
protected:
    template <typename T>
    VectorType _insert(const VectorType& a, const T& value, const int index) const noexcept
    {
        return _insert_helper(a, value, index, std::make_index_sequence<band_size> {});
    }
    VectorType _insert_bottom(VectorType a, const ScoreType value) const noexcept
    {
        return _insert<0>(a, value);
    }
    static VectorType _insert_top(VectorType a, const ScoreType value) noexcept
    {
        return _insert<num_blocks_ * block_words_ - 1>(a, value);
    }
private:
    static BlockType do_add(const BlockType& lhs, const BlockType& rhs, short) noexcept
    {
        return _mm_add_epi16(lhs, rhs);
    }
    static BlockType do_add(const BlockType& lhs, const BlockType& rhs, int) noexcept
    {
        return _mm_add_epi32(lhs, rhs);
    }
protected:
    static VectorType _add(const VectorType& lhs, const VectorType& rhs) noexcept
    {
        return transform([] (const auto& lhs, const auto& rhs) noexcept { return do_add(lhs, rhs, ScoreType {}); }, lhs, rhs);
    }
    static VectorType _and(const VectorType& lhs, const VectorType& rhs) noexcept
    {
        return transform([] (const auto& lhs, const auto& rhs) noexcept { return _mm_and_si128(lhs, rhs); }, lhs, rhs);
    }
    static VectorType _andnot(const VectorType& lhs, const VectorType& rhs) noexcept
    {
        return transform([] (const auto& lhs, const auto& rhs) noexcept { return _mm_andnot_si128(lhs, rhs); }, lhs, rhs);
    }
    static VectorType _or(const VectorType& lhs, const VectorType& rhs) noexcept
    {
        return transform([] (const auto& lhs, const auto& rhs) noexcept { return _mm_or_si128(lhs, rhs); }, lhs, rhs);
    }
private:
    static BlockType do_cmpeq(const BlockType& lhs, const BlockType& rhs, short) noexcept
    {
        return _mm_cmpeq_epi16(lhs, rhs);
    }
    static BlockType do_cmpeq(const BlockType& lhs, const BlockType& rhs, int) noexcept
    {
        return _mm_cmpeq_epi32(lhs, rhs);
    }
protected:
    static VectorType _cmpeq(const VectorType& lhs, const VectorType& rhs) noexcept
    {
        return transform([] (const auto& lhs, const auto& rhs) noexcept { return do_cmpeq(lhs, rhs, ScoreType {}); }, lhs, rhs);
    }
    static VectorType _left_shift_word(VectorType a) noexcept
    {
        adjacent_apply_reverse([] (const auto& lhs, const auto& rhs) noexcept {
            return _mm_or_si128(_mm_slli_si128(rhs, word_size), _mm_srli_si128(lhs, block_bytes_ - word_size));
        }, a, a);
        std::get<0>(a) = _mm_slli_si128(std::get<0>(a), word_size);
        return a;
    }
    static VectorType _right_shift_word(VectorType a) noexcept
    {
        adjacent_apply([] (const auto& lhs, const auto& rhs) noexcept {
            return _mm_or_si128(_mm_srli_si128(lhs, word_size), _mm_slli_si128(rhs, block_bytes_ - word_size));
        }, a, a);
        std::get<num_blocks - 1>(a) = _mm_srli_si128(std::get<num_blocks - 1>(a), word_size);
        return a;
    }
private:
    template <int n>
    static BlockType do_left_shift_bits(const BlockType& a, short) noexcept
    {
        return _mm_slli_epi16(a, n);
    }
    template <int n>
    static BlockType do_left_shift_bits(const BlockType& a, int) noexcept
    {
        return _mm_slli_epi32(a, n);
    }
protected:
    template <int n>
    static VectorType _left_shift_bits(const VectorType& a) noexcept
    {
        return transform([] (const auto& x) noexcept { return do_left_shift_bits<n>(x, ScoreType {}); }, a);
    }
private:
    template <int n>
    static BlockType do_right_shift_bits(const BlockType& a, short) noexcept
    {
        return _mm_srli_epi16(a, n);
    }
    template <int n>
    static BlockType do_right_shift_bits(const BlockType& a, int) noexcept
    {
        return _mm_srli_epi32(a, n);
    }
protected:
    template <int n>
    static VectorType _right_shift_bits(const VectorType& a) noexcept
    {
        return transform([] (const auto& x) noexcept { return do_right_shift_bits<n>(x, ScoreType {}); }, a);
    }
private:
    static BlockType do_min(const BlockType& lhs, const BlockType& rhs, short) noexcept
    {
        return _mm_min_epi16(lhs, rhs);
    }
    static BlockType do_min(const BlockType& lhs, const BlockType& rhs, int) noexcept
    {
        return _mm_min_epi32(lhs, rhs);
    }
protected:
    static VectorType _min(const VectorType& lhs, const VectorType& rhs) noexcept
    {
        return transform([] (const auto& lhs, const auto& rhs) noexcept { return do_min(lhs, rhs, ScoreType {}); }, lhs, rhs);
    }
private:
    static BlockType do_max(const BlockType& lhs, const BlockType& rhs, short) noexcept
    {
        return _mm_max_epi16(lhs, rhs);
    }
    static BlockType do_max(const BlockType& lhs, const BlockType& rhs, int) noexcept
    {
        return _mm_max_epi32(lhs, rhs);
    }
protected:
    static VectorType _max(const VectorType& lhs, const VectorType& rhs) noexcept
    {
        return transform([] (const auto& lhs, const auto& rhs) noexcept { return do_max(lhs, rhs, ScoreType {}); }, lhs, rhs);
    }
};

} // namespace simd
} // namespace hmm
} // namespace octopus

#endif
