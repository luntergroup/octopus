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

template <unsigned BandSize,
          typename ScoreTp = short>
class SSE2PairHMMInstructionSet
{
    using BlockType = __m128i;
    
protected:
    using ScoreType = ScoreTp;
    
    static_assert(is_short_or_int<ScoreType>, "ScoreType not short or int");
    
    constexpr static int word_size = sizeof(ScoreType);
    
    constexpr static const char* name = "SSE2";

private:
    constexpr static auto block_bytes_       = sizeof(BlockType);
    constexpr static auto block_words_       = block_bytes_ / word_size;
    constexpr static std::size_t num_blocks_ = BandSize / block_words_;
    
    static_assert(BandSize > 0, "BandSize must be positive");
    static_assert(BandSize % block_words_ == 0, "BandSize must be multiple of block words");
    
protected:
    using VectorType = std::conditional_t<num_blocks_ == 1, BlockType, std::array<BlockType, num_blocks_>>;
    
    constexpr static int band_size = num_blocks_ * block_words_;
    
    static_assert(sizeof(VectorType) / word_size == band_size, "size error");
    
private:
    static BlockType do_vectorise(ScoreType x, std::true_type, short) noexcept
    {
        return _mm_set1_epi16(x);
    }
    static BlockType do_vectorise(ScoreType x, std::true_type, int) noexcept
    {
        return _mm_set1_epi32(x);
    }
    static VectorType do_vectorise(ScoreType x, std::false_type, short) noexcept
    {
        return make_array<num_blocks_>(do_vectorise(x, std::true_type {}, short {}));
    }
    static VectorType do_vectorise(ScoreType x, std::false_type, int) noexcept
    {
        return make_array<num_blocks_>(do_vectorise(x, std::true_type {}, int {}));
    }
protected:
    static VectorType vectorise(ScoreType x) noexcept
    {
        return do_vectorise(x, std::is_same<VectorType, BlockType> {}, ScoreType {});
    }
private:
    template <typename T>
    static BlockType do_vectorise(const T* values, short) noexcept
    {
        return _mm_set_epi16(values[7], values[6], values[5], values[4],
                             values[3], values[2], values[1], values[0]);
    }
    template <typename T, std::size_t...Is>
    static BlockType do_vectorise(const T* values, int) noexcept
    {
        return _mm_set_epi32(values[3], values[2], values[1], values[0]);
    }
    template <typename T>
    static VectorType do_vectorise(const T* values, std::true_type) noexcept
    {
        return do_vectorise(values, ScoreType {});
    }
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
    template <typename T>
    static VectorType do_vectorise(const T* values, std::false_type) noexcept
    {
        return do_vectorise(values, std::make_index_sequence<num_blocks_>(), ScoreType {});
    }
protected:
    template <typename T>
    static VectorType vectorise(const T* values) noexcept
    {
        return do_vectorise(values, std::is_same<VectorType, BlockType> {});
    }
private:
    static BlockType do_vectorise_zero_set_last(ScoreType x, std::true_type, short) noexcept
    {
        return _mm_set_epi16(0,0,0,0,0,0,0,x);
    }
    static BlockType do_vectorise_zero_set_last(ScoreType x, std::true_type, int) noexcept
    {
        return _mm_set_epi32(0,0,0,x);
    }
    static VectorType do_vectorise_zero_set_last(ScoreType x, std::false_type, short) noexcept
    {
        return make_array<num_blocks_>(_mm_set_epi16(0,0,0,0,0,0,0,x), _mm_set_epi16(0,0,0,0,0,0,0,0));
    }
    static VectorType do_vectorise_zero_set_last(ScoreType x, std::false_type, int) noexcept
    {
        return make_array<num_blocks_>(_mm_set_epi32(0,0,0,x), _mm_set_epi32(0,0,0,0));
    }
protected:
    static VectorType vectorise_zero_set_last(ScoreType x) noexcept
    {
        return do_vectorise_zero_set_last(x, std::is_same<VectorType, BlockType> {},  ScoreType {});
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
    template <int index>
    static auto do_extract(const VectorType a, std::true_type) noexcept
    {
        return do_extract<index>(a, ScoreType {});
    }
    template <int index>
    static auto do_extract(const VectorType a, std::false_type) noexcept
    {
        constexpr static auto block_index = index / block_words_;
        static_assert(block_index < num_blocks_, "block index out range");
        constexpr static auto word_index = index % block_words_;
        static_assert(word_index < block_words_, "word index out range");
        return do_extract<word_index>(std::get<block_index>(a), ScoreType {});
    }
protected:
    template <int index>
    static auto _extract(const VectorType a) noexcept
    {
        return do_extract<index>(a, std::is_same<VectorType, BlockType> {});
    }
private:
    template <std::size_t... Is>
    static auto _extract_helper(const VectorType& a, const int index, std::index_sequence<Is...>) noexcept
    {
        if (index >= 0 && index < band_size) {
            decltype(_extract<0>(a)) result;
            int unused[] = {(index == Is ? (result = _extract<Is>(a), 0) : 0)...};
            (void) unused;
            return result;
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
    template <int index, typename T>
    static VectorType do_insert(const VectorType& a, T value, std::true_type) noexcept
    {
        return do_insert<index>(a, value, ScoreType {});
    }
    template <int index, typename T>
    static VectorType do_insert(VectorType a, T value, std::false_type) noexcept
    {
        constexpr static auto block_index = index / block_words_;
        static_assert(block_index < num_blocks_, "block index out range");
        constexpr static auto word_index = index % block_words_;
        static_assert(word_index < block_words_, "word index out range");
        std::get<block_index>(a) = do_insert<word_index>(std::get<block_index>(a), value, ScoreType {});
        return a;
    }
protected:
    template <int index, typename T>
    static VectorType _insert(const VectorType& a, const T& value) noexcept
    {
        return do_insert<index>(a, value, std::is_same<VectorType, BlockType> {});
    }
private:
    template <typename T, std::size_t... Is>
    static VectorType _insert_helper(const VectorType& a, const T& value, const int index, std::index_sequence<Is...>) noexcept
    {
        if (index >= 0 && index < band_size) {
            VectorType result;
            int unused[] = {(index == Is ? (result = _insert<Is>(a, value), 0) : 0)...};
            (void) unused;
            return result;
        } else {
            return a;
        }
    }
protected:
    template <typename T>
    static VectorType _insert(const VectorType& a, const T& value, const int index) noexcept
    {
        return _insert_helper(a, value, index, std::make_index_sequence<band_size> {});
    }
    static VectorType _insert_bottom(VectorType a, const ScoreType value) noexcept
    {
        return _insert<0>(a, value);
    }
    static VectorType _insert_top(VectorType a, const ScoreType value) noexcept
    {
        return _insert<band_size - 1>(a, value);
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
    static VectorType do_add(const VectorType& lhs, const VectorType& rhs, std::true_type) noexcept
    {
        return do_add(lhs, rhs, ScoreType {});
    }
    static VectorType do_add(const VectorType& lhs, const VectorType& rhs, std::false_type) noexcept
    {
        return transform([] (const auto& lhs, const auto& rhs) noexcept { return do_add(lhs, rhs, ScoreType {}); }, lhs, rhs);
    }
protected:
    static VectorType _add(const VectorType& lhs, const VectorType& rhs) noexcept
    {
        return do_add(lhs, rhs, std::is_same<VectorType, BlockType> {});
    }
private:
    static BlockType do_and(const BlockType& lhs, const BlockType& rhs, std::true_type) noexcept
    {
        return  _mm_and_si128(lhs, rhs);
    }
    static VectorType do_and(const VectorType& lhs, const VectorType& rhs, std::false_type) noexcept
    {
        return transform([] (const auto& lhs, const auto& rhs) noexcept { return do_and(lhs, rhs, std::true_type {}); }, lhs, rhs);
    }
    static BlockType do_andnot(const BlockType& lhs, const BlockType& rhs, std::true_type) noexcept
    {
        return  _mm_andnot_si128(lhs, rhs);
    }
    static VectorType do_andnot(const VectorType& lhs, const VectorType& rhs, std::false_type) noexcept
    {
        return transform([] (const auto& lhs, const auto& rhs) noexcept { return do_andnot(lhs, rhs, std::true_type {}); }, lhs, rhs);
    }
    static BlockType do_or(const BlockType& lhs, const BlockType& rhs, std::true_type) noexcept
    {
        return  _mm_or_si128(lhs, rhs);
    }
    static VectorType do_or(const VectorType& lhs, const VectorType& rhs, std::false_type) noexcept
    {
        return transform([] (const auto& lhs, const auto& rhs) noexcept { return do_or(lhs, rhs, std::true_type {}); }, lhs, rhs);
    }
protected:
    static VectorType _and(const VectorType& lhs, const VectorType& rhs) noexcept
    {
        return do_and(lhs, rhs, std::is_same<VectorType, BlockType> {});
    }
    static VectorType _andnot(const VectorType& lhs, const VectorType& rhs) noexcept
    {
        return do_andnot(lhs, rhs, std::is_same<VectorType, BlockType> {});
    }
    static VectorType _or(const VectorType& lhs, const VectorType& rhs) noexcept
    {
        return do_or(lhs, rhs, std::is_same<VectorType, BlockType> {});
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
    static VectorType do_cmpeq(const VectorType& lhs, const VectorType& rhs, std::true_type) noexcept
    {
        return do_cmpeq(lhs, rhs, ScoreType {});
    }
    static VectorType do_cmpeq(const VectorType& lhs, const VectorType& rhs, std::false_type) noexcept
    {
        return transform([] (const auto& lhs, const auto& rhs) noexcept { return do_cmpeq(lhs, rhs, ScoreType {}); }, lhs, rhs);
    }
protected:
    static VectorType _cmpeq(const VectorType& lhs, const VectorType& rhs) noexcept
    {
        return do_cmpeq(lhs, rhs, std::is_same<VectorType, BlockType> {});
    }
private:
    static VectorType do_left_shift_word(VectorType a, std::true_type) noexcept
    {
        return _mm_slli_si128(a, word_size);
    }
    static VectorType do_left_shift_word(VectorType a, std::false_type) noexcept
    {
        adjacent_apply_reverse([] (const auto& lhs, const auto& rhs) noexcept {
            return _mm_or_si128(_mm_slli_si128(rhs, word_size), _mm_srli_si128(lhs, block_bytes_ - word_size));
        }, a, a);
        std::get<0>(a) = _mm_slli_si128(std::get<0>(a), word_size);
        return a;
    }
    static VectorType do_right_shift_word(VectorType a, std::true_type) noexcept
    {
        return _mm_srli_si128(a, word_size);
    }
    static VectorType do_right_shift_word(VectorType a, std::false_type) noexcept
    {
        adjacent_apply([] (const auto& lhs, const auto& rhs) noexcept {
            return _mm_or_si128(_mm_srli_si128(lhs, word_size), _mm_slli_si128(rhs, block_bytes_ - word_size));
        }, a, a);
        std::get<num_blocks_ - 1>(a) = _mm_srli_si128(std::get<num_blocks_ - 1>(a), word_size);
        return a;
    }
protected:
    static VectorType _left_shift_word(VectorType a) noexcept
    {
        return do_left_shift_word(a, std::is_same<VectorType, BlockType> {});
    }
    static VectorType _right_shift_word(VectorType a) noexcept
    {
        return do_right_shift_word(a, std::is_same<VectorType, BlockType> {});
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
    template <int n>
    static VectorType do_left_shift_bits(const VectorType& a, std::true_type) noexcept
    {
        return do_left_shift_bits<n>(a, ScoreType {});
    }
    template <int n>
    static VectorType do_left_shift_bits(const VectorType& a, std::false_type) noexcept
    {
        return transform([] (const auto& x) noexcept { return do_left_shift_bits<n>(x, ScoreType {}); }, a);
    }
protected:
    template <int n>
    static VectorType _left_shift_bits(const VectorType& a) noexcept
    {
        return do_left_shift_bits<n>(a, std::is_same<VectorType, BlockType> {});
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
    template <int n>
    static VectorType do_right_shift_bits(const VectorType& a, std::true_type) noexcept
    {
        return do_right_shift_bits<n>(a, ScoreType {});
    }
    template <int n>
    static VectorType do_right_shift_bits(const VectorType& a, std::false_type) noexcept
    {
        return transform([] (const auto& x) noexcept { return do_right_shift_bits<n>(x, ScoreType {}); }, a);
    }
protected:
    template <int n>
    static VectorType _right_shift_bits(const VectorType& a) noexcept
    {
        return do_right_shift_bits<n>(a, std::is_same<VectorType, BlockType> {});
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
    static VectorType do_min(const VectorType& lhs, const VectorType& rhs, std::true_type) noexcept
    {
        return do_min(lhs, rhs, ScoreType {});
    }
    static VectorType do_min(const VectorType& lhs, const VectorType& rhs, std::false_type) noexcept
    {
        return transform([] (const auto& lhs, const auto& rhs) noexcept { return do_min(lhs, rhs, ScoreType {}); }, lhs, rhs);
    }
protected:
    static VectorType _min(const VectorType& lhs, const VectorType& rhs) noexcept
    {
        return do_min(lhs, rhs, std::is_same<VectorType, BlockType> {});
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
    static VectorType do_max(const VectorType& lhs, const VectorType& rhs, std::true_type) noexcept
    {
        return do_max(lhs, rhs, ScoreType {});
    }
    static VectorType do_max(const VectorType& lhs, const VectorType& rhs, std::false_type) noexcept
    {
        return transform([] (const auto& lhs, const auto& rhs) noexcept { return do_max(lhs, rhs, ScoreType {}); }, lhs, rhs);
    }
protected:
    static VectorType _max(const VectorType& lhs, const VectorType& rhs) noexcept
    {
        return do_max(lhs, rhs, std::is_same<VectorType, BlockType> {});
    }
};

} // namespace simd
} // namespace hmm
} // namespace octopus

#endif
