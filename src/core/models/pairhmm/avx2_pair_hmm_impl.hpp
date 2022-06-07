// Copyright (c) 2015-2021 Daniel Cooke and Gerton Lunter
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef avx2_pair_hmm_impl_hpp
#define avx2_pair_hmm_impl_hpp

#if __GNUC__ >= 6
#pragma GCC diagnostic ignored "-Wignored-attributes"
#endif

#include <cstdint>
#include <array>
#include <tuple>
#include <type_traits>
#include <cassert>
#if defined(__AVX2__)
#include <immintrin.h>
#endif // defined(__AVX2__)
#include "utils/array_tricks.hpp"
#include "system.hpp"

namespace octopus { namespace hmm { namespace simd {

#if defined(__AVX2__)

#define AVX2_PHMM

namespace detail {

// source: https://stackoverflow.com/a/25264853/2970186

template <int n>
std::enable_if_t<(n < 16), __m256i> _shift_left(const __m256i& a) noexcept
{
    return _mm256_alignr_epi8(a, _mm256_permute2x128_si256(a, a, _MM_SHUFFLE(0, 0, 2, 0)), 16 - n);
};
template <int n>
std::enable_if_t<(n == 16), __m256i> _shift_left(const __m256i& a) noexcept
{
    return _mm256_permute2x128_si256(a, a, _MM_SHUFFLE(0, 0, 2, 0));
}
template <int n>
std::enable_if_t<(n > 16), __m256i> _shift_left(const __m256i& a) noexcept
{
    return _shift_left<n - 16>(_mm256_permute2x128_si256(a, a, _MM_SHUFFLE(0, 0, 2, 0)));
};

template <int n>
std::enable_if_t<(n < 16), __m256i> _shift_right(const __m256i& a) noexcept
{
    return _mm256_alignr_epi8(_mm256_permute2x128_si256(a, a, _MM_SHUFFLE(2, 0, 0, 1)), a, n);
}
template <int n>
std::enable_if_t<(n == 16), __m256i> _shift_right(const __m256i& a) noexcept
{
    return _mm256_permute2x128_si256(a, a, _MM_SHUFFLE(2, 0, 0, 1));
}
template <int n>
std::enable_if_t<(n > 16), __m256i> _shift_right(const __m256i& a) noexcept
{
    return _shift_right<n - 16>(_mm256_permute2x128_si256(a, a, _MM_SHUFFLE(2, 0, 0, 1)));
}

} // namespace detail

template <unsigned BandSize = 16,
          typename ScoreTp = short>
class AVX2PairHMMInstructionSet
{
    using BlockType = __m256i;
    
protected:
    using ScoreType = ScoreTp;
    
    static_assert(std::is_same<ScoreType, short>::value || std::is_same<ScoreType, int>::value, "ScoreType not short or int");
    
    constexpr static int word_size = sizeof(ScoreType);
    
    constexpr static const char* name = "AVX2";

private:
    constexpr static auto block_bytes_       = sizeof(BlockType);
    constexpr static auto block_words_       = block_bytes_ / word_size;
    constexpr static std::size_t num_blocks_ = BandSize / block_words_;
    
    static_assert(BandSize > 0, "BandSize must be positive");
    static_assert(BandSize % block_words_ == 0, "BandSize must be multiple of block words");
    
protected:
    using VectorType = std::array<BlockType, num_blocks_>;
    
    constexpr static int band_size = num_blocks_ * block_words_;
    
    static_assert(sizeof(VectorType) / word_size == band_size, "size error");
    
private:
    static VectorType do_vectorise(ScoreType x, short) noexcept
    {
        return make_array<num_blocks_>(_mm256_set1_epi16(x));
    }
    static VectorType do_vectorise(ScoreType x, int) noexcept
    {
        return make_array<num_blocks_>(_mm256_set1_epi32(x));
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
            _mm256_set_epi16(values[block_words_ * Is + 15],
                             values[block_words_ * Is + 14],
                             values[block_words_ * Is + 13],
                             values[block_words_ * Is + 12],
                             values[block_words_ * Is + 11],
                             values[block_words_ * Is + 10],
                             values[block_words_ * Is + 9],
                             values[block_words_ * Is + 8],
                             values[block_words_ * Is + 7],
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
            _mm256_set_epi32(values[block_words_ * Is + 7],
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
    template <typename T>
    static VectorType vectorise(const T* values) noexcept
    {
        return do_vectorise(values, std::make_index_sequence<num_blocks_>(), ScoreType {});
    }
private:
    static VectorType do_vectorise_zero_set_last(ScoreType x, short) noexcept
    {
        return make_array<num_blocks_>(_mm256_set_epi16(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,x), _mm256_set_epi16(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0));
    }
    static VectorType do_vectorise_zero_set_last(ScoreType x, int) noexcept
    {
        return make_array<num_blocks_>(_mm256_set_epi32(0,0,0,0,0,0,0,x), _mm256_set_epi32(0,0,0,0,0,0,0,0));
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
        return _mm256_extract_epi16(a, index);
    }
    template <int index>
    static auto do_extract(const BlockType& a, int) noexcept
    {
        return _mm256_extract_epi32(a, index);
    }
protected:
    template <int index>
    static auto _extract(const VectorType a) noexcept
    {
        constexpr static auto block_index = index / block_words_;
        static_assert(block_index < num_blocks_, "block index out range");
        constexpr static auto word_index = index % block_words_;
        static_assert(word_index < block_words_, "word index out range");
        return do_extract<word_index>(std::get<block_index>(a), ScoreType {});
    }
private:
    template <std::size_t... Is>
    static auto _extract_helper(const VectorType& a, const int index, std::index_sequence<Is...>) noexcept
    {
        assert(index >= 0 && index < band_size);
        decltype(_extract<0>(a)) result;
        int unused[] = {(index == Is ? (result = _extract<Is>(a), 0) : 0)...};
        (void) unused;
        return result;
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
        return _mm256_insert_epi16(a, value, index);
    }
    template <int index, typename T>
    static BlockType do_insert(const BlockType& a, T value, int) noexcept
    {
        return _mm256_insert_epi32(a, value, index);
    }
protected:
    template <int index, typename T>
    static VectorType _insert(VectorType a, T value) noexcept
    {
        constexpr static auto block_index = index / block_words_;
        static_assert(block_index < num_blocks_, "block index out range");
        constexpr static auto word_index = index % block_words_;
        static_assert(word_index < block_words_, "word index out range");
        std::get<block_index>(a) = do_insert<word_index>(std::get<block_index>(a), value, ScoreType {});
        return a;
    }
private:
    template <typename T, std::size_t... Is>
    static VectorType _insert_helper(const VectorType& a, const T& value, const int index, std::index_sequence<Is...>) noexcept
    {
        assert(index >= 0);
        if (index < band_size) {
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
        return _mm256_add_epi16(lhs, rhs);
    }
    static BlockType do_add(const BlockType& lhs, const BlockType& rhs, int) noexcept
    {
        return _mm256_add_epi32(lhs, rhs);
    }
protected:
    static VectorType _add(const VectorType& lhs, const VectorType& rhs) noexcept
    {
        return transform([] (const auto& lhs, const auto& rhs) noexcept { return do_add(lhs, rhs, ScoreType {}); }, lhs, rhs);
    }
    static VectorType _and(const VectorType& lhs, const VectorType& rhs) noexcept
    {
        return transform([] (const auto& lhs, const auto& rhs) noexcept { return _mm256_and_si256(lhs, rhs); }, lhs, rhs);
    }
    static VectorType _andnot(const VectorType& lhs, const VectorType& rhs) noexcept
    {
        return transform([] (const auto& lhs, const auto& rhs) noexcept { return _mm256_andnot_si256(lhs, rhs); }, lhs, rhs);
    }
    static VectorType _or(const VectorType& lhs, const VectorType& rhs) noexcept
    {
        return transform([] (const auto& lhs, const auto& rhs) noexcept { return _mm256_or_si256(lhs, rhs); }, lhs, rhs);
    }
private:
    static BlockType do_cmpeq(const BlockType& lhs, const BlockType& rhs, short) noexcept
    {
        return _mm256_cmpeq_epi16(lhs, rhs);
    }
    static BlockType do_cmpeq(const BlockType& lhs, const BlockType& rhs, int) noexcept
    {
        return _mm256_cmpeq_epi32(lhs, rhs);
    }
protected:
    static VectorType _cmpeq(const VectorType& lhs, const VectorType& rhs) noexcept
    {
        return transform([] (const auto& lhs, const auto& rhs) noexcept { return do_cmpeq(lhs, rhs, ScoreType {}); }, lhs, rhs);
    }
    static VectorType _left_shift_word(VectorType a) noexcept
    {
        using detail::_shift_left;
        using detail::_shift_right;
        adjacent_apply_reverse([] (const auto& lhs, const auto& rhs) noexcept {
            return _mm256_or_si256(_shift_left<word_size>(rhs), _shift_right<block_bytes_ - word_size>(lhs));
        }, a, a);
        std::get<0>(a) = _shift_left<word_size>(std::get<0>(a));
        return a;
    }
    static VectorType _right_shift_word(VectorType a) noexcept
    {
        using detail::_shift_left;
        using detail::_shift_right;
        adjacent_apply([] (const auto& lhs, const auto& rhs) noexcept {
            return _mm256_or_si256(_shift_right<word_size>(lhs), _shift_left<block_bytes_ - word_size>(rhs));
        }, a, a);
        std::get<num_blocks_ - 1>(a) = _shift_right<word_size>(std::get<num_blocks_ - 1>(a));
        return a;
    }
private:
    template <int n>
    static BlockType do_left_shift_bits(const BlockType& a, short) noexcept
    {
        return _mm256_slli_epi16(a, n);
    }
    template <int n>
    static BlockType do_left_shift_bits(const BlockType& a, int) noexcept
    {
        return _mm256_slli_epi32(a, n);
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
        return _mm256_srli_epi16(a, n);
    }
    template <int n>
    static BlockType do_right_shift_bits(const BlockType& a, int) noexcept
    {
        return _mm256_srli_epi32(a, n);
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
        return _mm256_min_epi16(lhs, rhs);
    }
    static BlockType do_min(const BlockType& lhs, const BlockType& rhs, int) noexcept
    {
        return _mm256_min_epi32(lhs, rhs);
    }
protected:
    static VectorType _min(const VectorType& lhs, const VectorType& rhs) noexcept
    {
        return transform([] (const auto& lhs, const auto& rhs) noexcept { return do_min(lhs, rhs, ScoreType {}); }, lhs, rhs);
    }
private:
    static BlockType do_max(const BlockType& lhs, const BlockType& rhs, short) noexcept
    {
        return _mm256_max_epi16(lhs, rhs);
    }
    static BlockType do_max(const BlockType& lhs, const BlockType& rhs, int) noexcept
    {
        return _mm256_max_epi32(lhs, rhs);
    }
protected:
    static VectorType _max(const VectorType& lhs, const VectorType& rhs) noexcept
    {
        return transform([] (const auto& lhs, const auto& rhs) noexcept { return do_max(lhs, rhs, ScoreType {}); }, lhs, rhs);
    }
};

#endif // defined(__AVX2__)

} // namespace simd
} // namespace hmm
} // namespace octopus

#endif
