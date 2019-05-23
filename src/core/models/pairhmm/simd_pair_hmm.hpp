// Copyright (c) 2015-2019 Daniel Cooke and Gerton Lunter
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef simd_pair_hmm_hpp
#define simd_pair_hmm_hpp

#if __GNUC__ >= 6
    #pragma GCC diagnostic ignored "-Wignored-attributes"
#endif

#include <cstddef>
#include <algorithm>
#include <iterator>
#include <cassert>
#include <emmintrin.h>

#include <boost/container/small_vector.hpp>

namespace octopus { namespace hmm { namespace simd {

struct SSE2 {};
struct AVX2 {};
struct AVX512 {};

template <typename T>
struct ScoreType
{
    using type = T;
};
template <> struct ScoreType<SSE2>
{
    using type = short;
};
template <> struct ScoreType<AVX2>
{
    using type = short;
};

template <typename T>
using score_type_t = typename ScoreType<T>::type;

template <typename T>
struct VectorType
{
    using type = T;
};
template <> struct VectorType<SSE2>
{
    using type = __m128i;
};
template <> struct VectorType<AVX2>
{
    using type = std::int16_t;
};
// template <> struct VectorType<AVX512>
// {
//     using type = std::int32_t;
// };

template <typename T>
using vector_type_t = typename VectorType<T>::type;

inline vector_type_t<SSE2> vectorise(short penalty, SSE2) noexcept
{
    return _mm_set1_epi16(penalty);
}
// vector_type_t<AVX2> vectorise(short penalty, AVX2) noexcept
// {
//     return _mm_set1_epi16(penalty);
// }

inline vector_type_t<SSE2> vectorise_zero_set_last(short penalty, SSE2) noexcept
{
    return _mm_set_epi16(0,0,0,0,0,0,0,penalty);
}
inline vector_type_t<AVX2> vectorise_zero_set_last(short penalty, AVX2) noexcept
{
    return 0;
}

inline vector_type_t<SSE2> vectorise_reverse(const char* sequence, SSE2) noexcept
{
    return _mm_set_epi16(sequence[7], sequence[6], sequence[5], sequence[4],
                         sequence[3], sequence[2], sequence[1], sequence[0]);
}
// vector_type_t<AVX2> vectorise_reverse(const char* sequence, AVX2) noexcept
// {
//     return _mm_set1_epi16(penalty);
// }

inline vector_type_t<SSE2> vectorise_reverse_lshift(const std::int8_t* values, const int shift, SSE2) noexcept
{
    return _mm_set_epi16(values[7] << shift, values[6] << shift, values[5] << shift, values[4] << shift,
                         values[3] << shift, values[2] << shift, values[1] << shift, values[0] << shift);
}
inline vector_type_t<SSE2> vectorise_reverse_lshift(const std::int8_t value, const int shift, SSE2) noexcept
{
    return vectorise(value << shift, SSE2 {});
}

template <int imm>
inline auto simd_extract(const vector_type_t<SSE2> a) noexcept
{
    return _mm_extract_epi16(a, imm);
}

inline auto simd_extract(const vector_type_t<SSE2> a, const int n) noexcept
{
    switch (n) {
        case 0:  return simd_extract<0>(a);
        case 1:  return simd_extract<1>(a);
        case 2:  return simd_extract<2>(a);
        case 3:  return simd_extract<3>(a);
        case 4:  return simd_extract<4>(a);
        case 5:  return simd_extract<5>(a);
        case 6:  return simd_extract<6>(a);
        default: return simd_extract<7>(a);
    }
}

template <int imm, typename T>
inline vector_type_t<SSE2> simd_insert(const vector_type_t<SSE2>& a, T i) noexcept
{
    return _mm_insert_epi16(a, i, imm);
}
template <int imm, typename T>
inline vector_type_t<AVX2> simd_insert(const vector_type_t<AVX2>& a, T i) noexcept
{
    return a;
}

template <typename T>
inline vector_type_t<SSE2> simd_insert(const vector_type_t<SSE2>& a, const T i, const int n) noexcept
{
    switch (n) {
        case 0:  return simd_insert<0>(a, i);
        case 1:  return simd_insert<1>(a, i);
        case 2:  return simd_insert<2>(a, i);
        case 3:  return simd_insert<3>(a, i);
        case 4:  return simd_insert<4>(a, i);
        case 5:  return simd_insert<5>(a, i);
        case 6:  return simd_insert<6>(a, i);
        default: return simd_insert<7>(a, i);
    }
}

inline vector_type_t<SSE2> simd_add(const vector_type_t<SSE2>& lhs, const vector_type_t<SSE2>& rhs) noexcept
{
    return _mm_add_epi16(lhs, rhs);
}
inline vector_type_t<AVX2> simd_add(const vector_type_t<AVX2>& lhs, const vector_type_t<AVX2>& rhs) noexcept
{
    return lhs + rhs;
}

inline vector_type_t<SSE2> simd_and(const vector_type_t<SSE2>& lhs, const vector_type_t<SSE2>& rhs) noexcept
{
    return _mm_and_si128(lhs, rhs);
}
inline vector_type_t<AVX2> simd_and(const vector_type_t<AVX2>& lhs, const vector_type_t<AVX2>& rhs) noexcept
{
    return lhs & rhs;
}

inline vector_type_t<SSE2> simd_andnot(const vector_type_t<SSE2>& lhs, const vector_type_t<SSE2>& rhs) noexcept
{
    return _mm_andnot_si128(lhs, rhs);
}
inline vector_type_t<AVX2> simd_andnot(const vector_type_t<AVX2>& lhs, const vector_type_t<AVX2>& rhs) noexcept
{
    return lhs & !rhs;
}

inline vector_type_t<SSE2> simd_or(const vector_type_t<SSE2>& lhs, const vector_type_t<SSE2>& rhs) noexcept
{
    return _mm_or_si128(lhs, rhs);
}
inline vector_type_t<AVX2> simd_or(const vector_type_t<AVX2>& lhs, const vector_type_t<AVX2>& rhs) noexcept
{
    return lhs | rhs;
}

inline vector_type_t<SSE2> simd_cmpeq(const vector_type_t<SSE2>& lhs, const vector_type_t<SSE2>& rhs) noexcept
{
    return _mm_cmpeq_epi16(lhs, rhs);
}
inline vector_type_t<AVX2> simd_cmpeq(const vector_type_t<AVX2>& lhs, const vector_type_t<AVX2>& rhs) noexcept
{
    return lhs == rhs;
}

template <int imm>
inline vector_type_t<SSE2> simd_slli_si(const vector_type_t<SSE2>& a) noexcept
{
    return _mm_slli_si128(a, imm);
}
template <int imm>
inline vector_type_t<AVX2> simd_slli_si(const vector_type_t<AVX2>& a) noexcept
{
    return a;
}
template <int imm>
inline vector_type_t<SSE2> simd_slli_epi(const vector_type_t<SSE2>& a) noexcept
{
    return _mm_slli_epi16(a, imm);
}
template <int imm>
inline vector_type_t<AVX2> simd_slli_epi(const vector_type_t<AVX2>& a) noexcept
{
    return a;
}
template <int imm>
inline vector_type_t<SSE2> simd_srli_si(const vector_type_t<SSE2>& a) noexcept
{
    return _mm_srli_si128(a, imm);
}
template <int imm>
inline vector_type_t<AVX2> simd_srli_si(const vector_type_t<AVX2>& a) noexcept
{
    return a;
}

inline vector_type_t<SSE2> simd_min(const vector_type_t<SSE2>& lhs, const vector_type_t<SSE2>& rhs) noexcept
{
    return _mm_min_epi16(lhs, rhs);
}
inline vector_type_t<AVX2> simd_min(const vector_type_t<AVX2>& lhs, const vector_type_t<AVX2>& rhs) noexcept
{
    return std::min(lhs, rhs);
}

inline vector_type_t<SSE2> simd_max(const vector_type_t<SSE2>& lhs, const vector_type_t<SSE2>& rhs) noexcept
{
    return _mm_max_epi16(lhs, rhs);
}
inline vector_type_t<AVX2> simd_max(const vector_type_t<AVX2>& lhs, const vector_type_t<AVX2>& rhs) noexcept
{
    return std::max(lhs, rhs);
}

constexpr int bit_shift(SSE2) { return 2; }
constexpr int bit_shift(AVX2) { return 2; }

constexpr int band_size(SSE2) { return 8; }
constexpr int band_size(AVX2) { return 16; }

constexpr score_type_t<SSE2> n_score(SSE2) { return 2 << bit_shift(SSE2 {}); }
constexpr score_type_t<AVX2> n_score(AVX2) { return 2 << bit_shift(AVX2 {}); }

constexpr score_type_t<SSE2> infinity(SSE2) { return 0x7800; }
constexpr score_type_t<AVX2> infinity(AVX2) { return 0x7800; }

template <typename SIMD>
inline void update_gap_penalty(vector_type_t<SIMD>& current, const std::int8_t* source, const std::size_t gap_idx) noexcept
{
    constexpr SIMD simd {};
    current = simd_insert(simd_srli_si<bit_shift(simd)>(current), source[gap_idx] << bit_shift(simd), band_size(simd) - 1);
}

template <typename SIMD>
inline void update_gap_penalty(vector_type_t<SIMD>& current, const short source, const std::size_t gap_idx) noexcept {}

template <typename SIMD,
          typename OpenPenalty,
          typename ExtendPenalty>
int
align(const char* truth,
      const char* target,
      const std::int8_t* qualities,
      const int truth_len,
      const int target_len,
      const OpenPenalty gap_open,
      const ExtendPenalty gap_extend,
      short nuc_prior) noexcept
{
    constexpr SIMD simd {};
    constexpr static auto inf = infinity(simd);
    constexpr static auto nscore = n_score(simd);
    constexpr static auto bandsize = band_size(simd);
    constexpr static auto lshift = bit_shift(simd);
    constexpr static int qual_factor {64}; // what does 64 mean?
    constexpr static int bitmask {0x8000}; // what does 64 mean?
    const static auto _inf = vectorise(inf, simd);
    const static auto _nscore_m_inf = vectorise(nscore - inf, simd);
    const static auto _n = vectorise('N', simd);
    assert(truth_len > band_size(simd) && (truth_len == target_len + 2 * band_size(simd) - 1));
    
    auto _m1 = vectorise(inf, simd);
    auto _i1 = _m1, _d1 = _m1, _m2 = _m1, _i2 = _m1, _d2 = _m1;
    const auto _nuc_prior  = vectorise(nuc_prior << lshift, simd);
    auto _initmask  = vectorise_zero_set_last(-1, simd);
    auto _initmask2 = vectorise_zero_set_last(-bitmask, simd);
    // truth is initialized with the n-long prefix, in forward direction
    // target is initialized as empty; reverse direction
    auto _truthwin     = vectorise_reverse(truth, simd);
    auto _targetwin    = _m1;
    auto _qualitieswin = vectorise(qual_factor << lshift, simd);
    auto _gap_open     = vectorise_reverse_lshift(gap_open, lshift, simd);
    auto _gap_extend   = vectorise_reverse_lshift(gap_extend, lshift, simd);
    auto _truthnqual   = simd_add(simd_and(simd_cmpeq(_truthwin, _n), _nscore_m_inf), _inf); // if N, make nScore; if != N, make inf
    
    score_type_t<SIMD> minscore {inf};
    for (int s {0}; s <= 2 * (target_len + bandsize); s += 2) {
        // truth is current; target needs updating
        _targetwin    = simd_slli_si<lshift>(_targetwin);
        _qualitieswin = simd_slli_si<lshift>(_qualitieswin);
        if (s / 2 < target_len) {
            _targetwin    = simd_insert(_targetwin, target[s / 2], 0);
            _qualitieswin = simd_insert(_qualitieswin, qualities[s / 2] << lshift, 0);
        } else {
            _targetwin    = simd_insert(_targetwin, '0', 0);
            _qualitieswin = simd_insert(_qualitieswin, qual_factor << lshift, 0);
        }
        // S even
        _m1 = simd_or(_initmask2, simd_andnot(_initmask, _m1));
        _m2 = simd_or(_initmask2, simd_andnot(_initmask, _m2));
        _m1 = simd_min(_m1, simd_min(_i1, _d1));
        if (s / 2 >= target_len) {
            minscore = std::min(static_cast<decltype(minscore)>(simd_extract(_m1, std::max(0, s / 2 - target_len))), minscore);
        }
        _m1 = simd_add(_m1, simd_min(simd_andnot(simd_cmpeq(_targetwin, _truthwin), _qualitieswin), _truthnqual));
        _d1 = simd_min(simd_add(_d2, _gap_extend), simd_add(simd_min(_m2, _i2), simd_srli_si<lshift>(_gap_open))); // allow I->D
        _d1 = simd_insert(simd_slli_si<lshift>(_d1), inf, 0);
        _i1 = simd_add(simd_min(simd_add(_i2, _gap_extend), simd_add(_m2, _gap_open)), _nuc_prior);
        // S odd
        // truth needs updating; target is current
        const auto pos = bandsize + s / 2;
        const bool pos_in_range {pos < truth_len};
        const char base {pos_in_range ? truth[pos] : 'N'};
        _truthwin   = simd_insert(simd_srli_si<lshift>(_truthwin), base, bandsize - 1);
        _truthnqual = simd_insert(simd_srli_si<lshift>(_truthnqual), base == 'N' ? nscore : inf, bandsize - 1);
        const auto gap_idx = pos_in_range ? pos : truth_len - 1;
        update_gap_penalty<SIMD>(_gap_open, gap_open, gap_idx);
        update_gap_penalty<SIMD>(_gap_extend, gap_extend, gap_idx);
        _initmask  = simd_slli_si<lshift>(_initmask);
        _initmask2 = simd_slli_si<lshift>(_initmask2);
        _m2 = simd_min(_m2, simd_min(_i2, _d2));
        if (s / 2 >= target_len) {
            minscore = std::min(static_cast<decltype(minscore)>(simd_extract(_m2, s / 2 - target_len)), minscore);
        }
        _m2 = simd_add(_m2, simd_min(simd_andnot(simd_cmpeq(_targetwin, _truthwin), _qualitieswin), _truthnqual));
        _d2 = simd_min(simd_add(_d1, _gap_extend), simd_add(simd_min(_m1, _i1), _gap_open)); // allow I->D
        _i2 = simd_insert(simd_add(simd_min(simd_add(simd_srli_si<lshift>(_i1), _gap_extend), simd_add(simd_srli_si<lshift>(_m1), _gap_open)), _nuc_prior), inf, bandsize - 1);
    }
    return (minscore + bitmask) >> lshift;
}

constexpr std::size_t staticBackpointerCapacity {10000};

template <typename T>
using SmallVector = boost::container::small_vector<T, staticBackpointerCapacity>;

constexpr static char gap_label {'-'};

template <typename SIMD,
          typename OpenPenalty,
          typename ExtendPenalty>
int
align(const char* truth,
      const char* target,
      const std::int8_t* qualities,
      const int truth_len,
      const int target_len,
      const OpenPenalty gap_open,
      const ExtendPenalty gap_extend,
      short nuc_prior,
      int& first_pos,
      char* align1,
      char* align2) noexcept
{
    constexpr SIMD simd {};
    constexpr static auto inf = infinity(simd);
    constexpr static auto nscore = n_score(simd);
    constexpr static auto bandsize = band_size(simd);
    constexpr static auto lshift = bit_shift(simd);
    constexpr static int qual_factor {64}; // what does 64 mean?
    constexpr static int bitmask {0x8000};
    const static auto _inf = vectorise(inf, simd);
    const static auto _nscore_m_inf = vectorise(nscore - inf, simd);
    const static auto _n = vectorise('N', simd);
    constexpr static int match_label  {0};
    constexpr static int insert_label {1};
    constexpr static int delete_label {3};
    static const auto _three = vectorise(3, simd);
    assert(truth_len > band_size(simd) && (truth_len == target_len + 2 * band_size(simd) - 1));
    
    auto _m1 = vectorise(inf, simd);
    auto _i1 = _m1, _d1 = _m1, _m2 = _m1, _i2 = _m1, _d2 = _m1;
    const auto _nuc_prior = vectorise(nuc_prior << lshift, simd);
    auto _initmask  = vectorise_zero_set_last(-1, simd);
    auto _initmask2 = vectorise_zero_set_last(-bitmask, simd);
    // truth is initialized with the n-long prefix, in forward direction
    // target is initialized as empty; reverse direction
    auto _truthwin     = vectorise_reverse(truth, simd);
    auto _targetwin    = _m1;
    auto _qualitieswin = vectorise(qual_factor << lshift, simd);
    auto _gap_open     = vectorise_reverse_lshift(gap_open, lshift, simd);
    auto _gap_extend   = vectorise_reverse_lshift(gap_extend, lshift, simd);
    auto _truthnqual   = simd_add(simd_and(simd_cmpeq(_truthwin, _n), _nscore_m_inf), _inf); // if N, make nScore; if != N, make inf
    
    SmallVector<vector_type_t<SIMD>> _backpointers(2 * (truth_len + bandsize));
    
    score_type_t<SIMD> minscore {inf}, cur_score;
    int s, minscoreidx {-1};
    for (s = 0; s <= 2 * (target_len + bandsize); s += 2) {
        // truth is current; target needs updating
        _targetwin    = simd_slli_si<lshift>(_targetwin);
        _qualitieswin = simd_slli_si<lshift>(_qualitieswin);
        if (s / 2 < target_len) {
            _targetwin    = simd_insert(_targetwin, target[s / 2], 0);
            _qualitieswin = simd_insert(_qualitieswin, qualities[s / 2] << lshift, 0);
        } else {
            _targetwin    = simd_insert(_targetwin, '0', 0);
            _qualitieswin = simd_insert(_qualitieswin, qual_factor << lshift, 0);
        }
        // S even
        _m1 = simd_or(_initmask2, simd_andnot(_initmask, _m1));
        _m2 = simd_or(_initmask2, simd_andnot(_initmask, _m2));
        _m1 = simd_min(_m1, simd_min(_i1, _d1));
        if (s / 2 >= target_len) {
            cur_score = simd_extract(_m1, s / 2 - target_len);
            if (cur_score < minscore) {
                minscore = cur_score;
                minscoreidx = s;
            }
        }
        _m1 = simd_add(_m1, simd_min(simd_andnot(simd_cmpeq(_targetwin, _truthwin), _qualitieswin), _truthnqual));
        _d1 = simd_min(simd_add(_d2, _gap_extend), simd_add(simd_min(_m2, _i2), simd_srli_si<lshift>(_gap_open))); // allow I->D
        _d1 = simd_insert(simd_slli_si<lshift>(_d1), inf, 0);
        _i1 = simd_add(simd_min(simd_add(_i2, _gap_extend), simd_add(_m2, _gap_open)), _nuc_prior);
        _backpointers[s] = simd_or(simd_or(simd_and(_three, _m1), simd_slli_epi<2 * insert_label>(simd_and(_three, _i1))),
                                   simd_slli_epi<2 * delete_label>(simd_and(_three, _d1)));
        // set state labels
        _m1 = simd_andnot(_three, _m1);
        _i1 = simd_or(simd_andnot(_three, _i1), simd_srli_si<1>(_three));
        _d1 = simd_or(simd_andnot(_three, _d1), _three);
        // S odd
        // truth needs updating; target is current
        const auto pos = bandsize + s / 2;
        const bool pos_in_range {pos < truth_len};
        const char base {pos_in_range ? truth[pos] : 'N'};
        _truthwin   = simd_insert(simd_srli_si<lshift>(_truthwin), base, bandsize - 1);
        _truthnqual = simd_insert(simd_srli_si<lshift>(_truthnqual), base == 'N' ? nscore : inf, bandsize - 1);
        const auto gap_idx = pos_in_range ? pos : truth_len - 1;
        update_gap_penalty<SIMD>(_gap_open, gap_open, gap_idx);
        update_gap_penalty<SIMD>(_gap_extend, gap_extend, gap_idx);
        _initmask  = simd_slli_si<lshift>(_initmask);
        _initmask2 = simd_slli_si<lshift>(_initmask2);
        _m2 = simd_min(_m2, simd_min(_i2, _d2));
        if (s / 2 >= target_len) {
            cur_score = simd_extract(_m2, s / 2 - target_len);
            if (cur_score < minscore) {
                minscore = cur_score;
                minscoreidx = s + 1;
            }
        }
        _m2 = simd_add(_m2, simd_min(simd_andnot(simd_cmpeq(_targetwin, _truthwin), _qualitieswin), _truthnqual));
        _d2 = simd_min(simd_add(_d1, _gap_extend), simd_add(simd_min(_m1, _i1), _gap_open)); // allow I->D
        _i2 = simd_insert(simd_add(simd_min(simd_add(simd_srli_si<lshift>(_i1), _gap_extend), simd_add(simd_srli_si<lshift>(_m1), _gap_open)), _nuc_prior), inf, bandsize - 1);
        _backpointers[s + 1] = simd_or(simd_or(simd_and(_three, _m2), simd_slli_epi<2 * insert_label>(simd_and(_three, _i2))),
                                       simd_slli_epi<2 * delete_label>(simd_and(_three, _d2)));
        // set state labels
        _m2 = simd_andnot(_three, _m2);
        _i2 = simd_or(simd_andnot(_three, _i2), simd_srli_si<1>(_three));
        _d2 = simd_or(simd_andnot(_three, _d2), _three);
    }
    if (minscoreidx < 0) {
        // minscore was never updated so we must have overflowed badly
        first_pos = -1;
        return -1;
    }
    s = minscoreidx;    // point to the dummy match transition
    auto i      = s / 2 - target_len;
    auto y      = target_len;
    auto x      = s - y;
    auto alnidx = 0;
    auto ptr = reinterpret_cast<short*>(_backpointers.data() + s);
    if ((ptr + i) < reinterpret_cast<short*>(_backpointers.data())
        || (ptr + i) >= reinterpret_cast<short*>(_backpointers.data() + _backpointers.size())) {
        first_pos = -1;
        return -1;
    }
    auto state = (ptr[i] >> (2 * match_label)) & 3;
    s -= 2;
    // this is 2*y (s even) or 2*y+1 (s odd)
    while (y > 0) {
        if (s < 0 || i < 0) {
            // This should never happen so must have overflowed
            first_pos = -1;
            return -1;
        }
        const auto new_state = (reinterpret_cast<short*>(_backpointers.data() + s)[i] >> (2 * state)) & 3;
        if (state == match_label) {
            s -= 2;
            align1[alnidx] = truth[--x];
            align2[alnidx] = target[--y];
        } else if (state == insert_label) {
            i += s & 1;
            s -= 1;
            align1[alnidx] = gap_label;
            align2[alnidx] = target[--y];
        } else {
            s -= 1;
            i -= s & 1;
            align1[alnidx] = truth[--x];
            align2[alnidx] = gap_label;
        }
        state = new_state;
        alnidx++;
    }
    align1[alnidx] = 0;
    align2[alnidx] = 0;
    first_pos = x;
    // reverse them
    for (int j {alnidx - 1}, i = 0; i < j; ++i, j--) {
        x = align1[i];
        y = align2[i];
        align1[i] = align1[j];
        align2[i] = align2[j];
        align1[j] = x;
        align2[j] = y;
    }
    return (minscore + bitmask) >> lshift;
}

} // namespace simd
} // namespace hmm
} // namespace octopus

#endif
