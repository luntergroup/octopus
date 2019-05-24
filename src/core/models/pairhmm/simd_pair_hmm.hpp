// Copyright (c) 2015-2019 Daniel Cooke and Gerton Lunter
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef _pair_hmm_hpp
#define _pair_hmm_hpp

#if __GNUC__ >= 6
    #pragma GCC diagnostic ignored "-Wignored-attributes"
#endif

#include <cstddef>
#include <algorithm>
#include <iterator>
#include <cassert>
#include <emmintrin.h>
#include <immintrin.h>

#include <boost/container/small_vector.hpp>

namespace octopus { namespace hmm { namespace simd {

template <typename InstructionSetPolicy>
class PairHMM : private InstructionSetPolicy
{
    // Types
    using VectorType = typename InstructionSetPolicy::VectorType;
    using ScoreType  = typename InstructionSetPolicy::ScoreType;
    using SmallVector = typename InstructionSetPolicy::SmallVector;
    // Constants
    using InstructionSetPolicy::band_size;
    using InstructionSetPolicy::infinity;
    using InstructionSetPolicy::trace_bits;
    using InstructionSetPolicy::n_score;
    // Methods
    using InstructionSetPolicy::vectorise;
    using InstructionSetPolicy::vectorise_reverse;
    using InstructionSetPolicy::vectorise_reverse_lshift;
    using InstructionSetPolicy::vectorise_zero_set_last;
    using InstructionSetPolicy::_extract;
    using InstructionSetPolicy::_insert;
    using InstructionSetPolicy::_add;
    using InstructionSetPolicy::_and;
    using InstructionSetPolicy::_andnot;
    using InstructionSetPolicy::_or;
    using InstructionSetPolicy::_cmpeq;
    using InstructionSetPolicy::_min;
    using InstructionSetPolicy::_max;
    
    template <int idx>
    auto _slli_si(const VectorType& vec) const noexcept { return InstructionSetPolicy::template _slli_si<idx>(vec); }
    template <int idx>
    auto _srli_si(const VectorType& vec) const noexcept { return InstructionSetPolicy::template _srli_si<idx>(vec); }
    template <int idx>
    auto _slli_epi(const VectorType& vec) const noexcept { return InstructionSetPolicy::template _slli_epi<idx>(vec); }
    
    void update_gap_penalty(VectorType& current, const std::int8_t* source, const std::size_t gap_idx) const noexcept
    {
        current = _insert(_srli_si<trace_bits>(current), source[gap_idx] << trace_bits, band_size - 1);
    }
    void update_gap_penalty(VectorType& current, const short source, const std::size_t gap_idx) const noexcept {}

public:
    constexpr static char gap_label = '-';
    
    template <typename OpenPenalty,
              typename ExtendPenalty>
    int align(const char* truth,
              const char* target,
              const std::int8_t* qualities,
              const int truth_len,
              const int target_len,
              const OpenPenalty gap_open,
              const ExtendPenalty gap_extend,
              short nuc_prior) const noexcept
    {
        constexpr static int qual_factor {64}; // what does 64 mean?
        constexpr static int bitmask {0x8000}; // ?
        const static auto _inf = vectorise(infinity);
        const static auto _nscore_m_inf = vectorise(n_score - infinity);
        const static auto _n = vectorise('N');
        assert(truth_len > band_size && (truth_len == target_len + 2 * band_size - 1));
        
        auto _m1 = vectorise(infinity);
        auto _i1 = _m1, _d1 = _m1, _m2 = _m1, _i2 = _m1, _d2 = _m1;
        const auto _nuc_prior  = vectorise(nuc_prior << trace_bits);
        auto _initmask  = vectorise_zero_set_last(-1);
        auto _initmask2 = vectorise_zero_set_last(-bitmask);
        // truth is initialized with the n-long prefix, in forward direction
        // target is initialized as empty; reverse direction
        auto _truthwin     = vectorise_reverse(truth);
        auto _targetwin    = _m1;
        auto _qualitieswin = vectorise(qual_factor << trace_bits);
        auto _gap_open     = vectorise_reverse_lshift(gap_open, trace_bits);
        auto _gap_extend   = vectorise_reverse_lshift(gap_extend, trace_bits);
        auto _truthnqual   = _add(_and(_cmpeq(_truthwin, _n), _nscore_m_inf), _inf); // if N, make nScore; if != N, make inf
        
        ScoreType minscore {infinity};
        for (int s {0}; s <= 2 * (target_len + band_size); s += 2) {
            // truth is current; target needs updating
            _targetwin    = _slli_si<trace_bits>(_targetwin);
            _qualitieswin = _slli_si<trace_bits>(_qualitieswin);
            if (s / 2 < target_len) {
                _targetwin    = _insert(_targetwin, target[s / 2], 0);
                _qualitieswin = _insert(_qualitieswin, qualities[s / 2] << trace_bits, 0);
            } else {
                _targetwin    = _insert(_targetwin, '0', 0);
                _qualitieswin = _insert(_qualitieswin, qual_factor << trace_bits, 0);
            }
            // S even
            _m1 = _or(_initmask2, _andnot(_initmask, _m1));
            _m2 = _or(_initmask2, _andnot(_initmask, _m2));
            _m1 = _min(_m1, _min(_i1, _d1));
            if (s / 2 >= target_len) {
                minscore = std::min(static_cast<decltype(minscore)>(_extract(_m1, std::max(0, s / 2 - target_len))), minscore);
            }
            _m1 = _add(_m1, _min(_andnot(_cmpeq(_targetwin, _truthwin), _qualitieswin), _truthnqual));
            _d1 = _min(_add(_d2, _gap_extend), _add(_min(_m2, _i2), _srli_si<trace_bits>(_gap_open))); // allow I->D
            _d1 = _insert(_slli_si<trace_bits>(_d1), infinity, 0);
            _i1 = _add(_min(_add(_i2, _gap_extend), _add(_m2, _gap_open)), _nuc_prior);
            // S odd
            // truth needs updating; target is current
            const auto pos = band_size + s / 2;
            const bool pos_in_range {pos < truth_len};
            const char base {pos_in_range ? truth[pos] : 'N'};
            _truthwin   = _insert(_srli_si<trace_bits>(_truthwin), base, band_size - 1);
            _truthnqual = _insert(_srli_si<trace_bits>(_truthnqual), base == 'N' ? n_score : infinity, band_size - 1);
            const auto gap_idx = pos_in_range ? pos : truth_len - 1;
            update_gap_penalty(_gap_open, gap_open, gap_idx);
            update_gap_penalty(_gap_extend, gap_extend, gap_idx);
            _initmask  = _slli_si<trace_bits>(_initmask);
            _initmask2 = _slli_si<trace_bits>(_initmask2);
            _m2 = _min(_m2, _min(_i2, _d2));
            if (s / 2 >= target_len) {
                minscore = std::min(static_cast<decltype(minscore)>(_extract(_m2, s / 2 - target_len)), minscore);
            }
            _m2 = _add(_m2, _min(_andnot(_cmpeq(_targetwin, _truthwin), _qualitieswin), _truthnqual));
            _d2 = _min(_add(_d1, _gap_extend), _add(_min(_m1, _i1), _gap_open)); // allow I->D
            _i2 = _insert(_add(_min(_add(_srli_si<trace_bits>(_i1), _gap_extend), _add(_srli_si<trace_bits>(_m1), _gap_open)), _nuc_prior), infinity, band_size - 1);
        }
        return (minscore + bitmask) >> trace_bits;
    }
    
    template <typename OpenPenalty,
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
        constexpr static int qual_factor {64}; // what does 64 mean?
        constexpr static int bitmask {0x8000}; // what does 64 mean?
        const static auto _inf = vectorise(infinity);
        const static auto _nscore_m_inf = vectorise(n_score - infinity);
        const static auto _n = vectorise('N');
        constexpr static int match_label  {0};
        constexpr static int insert_label {1};
        constexpr static int delete_label {3};
        static const auto _three = vectorise(3);
        assert(truth_len > band_size && (truth_len == target_len + 2 * band_size - 1));
        
        auto _m1 = vectorise(infinity);
        auto _i1 = _m1, _d1 = _m1, _m2 = _m1, _i2 = _m1, _d2 = _m1;
        const auto _nuc_prior = vectorise(nuc_prior << trace_bits);
        auto _initmask  = vectorise_zero_set_last(-1);
        auto _initmask2 = vectorise_zero_set_last(-bitmask);
        // truth is initialized with the n-long prefix, in forward direction
        // target is initialized as empty; reverse direction
        auto _truthwin     = vectorise_reverse(truth);
        auto _targetwin    = _m1;
        auto _qualitieswin = vectorise(qual_factor << trace_bits);
        auto _gap_open     = vectorise_reverse_lshift(gap_open, trace_bits);
        auto _gap_extend   = vectorise_reverse_lshift(gap_extend, trace_bits);
        auto _truthnqual   = _add(_and(_cmpeq(_truthwin, _n), _nscore_m_inf), _inf); // if N, make nScore; if != N, make inf
        
        SmallVector _backpointers(2 * (truth_len + band_size));
        
        ScoreType minscore {infinity}, cur_score;
        int s, minscoreidx {-1};
        for (s = 0; s <= 2 * (target_len + band_size); s += 2) {
            // truth is current; target needs updating
            _targetwin    = _slli_si<trace_bits>(_targetwin);
            _qualitieswin = _slli_si<trace_bits>(_qualitieswin);
            if (s / 2 < target_len) {
                _targetwin    = _insert(_targetwin, target[s / 2], 0);
                _qualitieswin = _insert(_qualitieswin, qualities[s / 2] << trace_bits, 0);
            } else {
                _targetwin    = _insert(_targetwin, '0', 0);
                _qualitieswin = _insert(_qualitieswin, qual_factor << trace_bits, 0);
            }
            // S even
            _m1 = _or(_initmask2, _andnot(_initmask, _m1));
            _m2 = _or(_initmask2, _andnot(_initmask, _m2));
            _m1 = _min(_m1, _min(_i1, _d1));
            if (s / 2 >= target_len) {
                cur_score = _extract(_m1, s / 2 - target_len);
                if (cur_score < minscore) {
                    minscore = cur_score;
                    minscoreidx = s;
                }
            }
            _m1 = _add(_m1, _min(_andnot(_cmpeq(_targetwin, _truthwin), _qualitieswin), _truthnqual));
            _d1 = _min(_add(_d2, _gap_extend), _add(_min(_m2, _i2), _srli_si<trace_bits>(_gap_open))); // allow I->D
            _d1 = _insert(_slli_si<trace_bits>(_d1), infinity, 0);
            _i1 = _add(_min(_add(_i2, _gap_extend), _add(_m2, _gap_open)), _nuc_prior);
            _backpointers[s] = _or(_or(_and(_three, _m1), _slli_epi<2 * insert_label>(_and(_three, _i1))),
                                            _slli_epi<2 * delete_label>(_and(_three, _d1)));
            // set state labels
            _m1 = _andnot(_three, _m1);
            _i1 = _or(_andnot(_three, _i1), _srli_si<1>(_three));
            _d1 = _or(_andnot(_three, _d1), _three);
            // S odd
            // truth needs updating; target is current
            const auto pos = band_size + s / 2;
            const bool pos_in_range {pos < truth_len};
            const char base {pos_in_range ? truth[pos] : 'N'};
            _truthwin   = _insert(_srli_si<trace_bits>(_truthwin), base, band_size - 1);
            _truthnqual = _insert(_srli_si<trace_bits>(_truthnqual), base == 'N' ? n_score : infinity, band_size - 1);
            const auto gap_idx = pos_in_range ? pos : truth_len - 1;
            update_gap_penalty(_gap_open, gap_open, gap_idx);
            update_gap_penalty(_gap_extend, gap_extend, gap_idx);
            _initmask  = _slli_si<trace_bits>(_initmask);
            _initmask2 = _slli_si<trace_bits>(_initmask2);
            _m2 = _min(_m2, _min(_i2, _d2));
            if (s / 2 >= target_len) {
                cur_score = _extract(_m2, s / 2 - target_len);
                if (cur_score < minscore) {
                    minscore = cur_score;
                    minscoreidx = s + 1;
                }
            }
            _m2 = _add(_m2, _min(_andnot(_cmpeq(_targetwin, _truthwin), _qualitieswin), _truthnqual));
            _d2 = _min(_add(_d1, _gap_extend), _add(_min(_m1, _i1), _gap_open)); // allow I->D
            _i2 = _insert(_add(_min(_add(_srli_si<trace_bits>(_i1), _gap_extend), _add(_srli_si<trace_bits>(_m1), _gap_open)), _nuc_prior), infinity, band_size - 1);
            _backpointers[s + 1] = _or(_or(_and(_three, _m2), _slli_epi<2 * insert_label>(_and(_three, _i2))),
                                                _slli_epi<2 * delete_label>(_and(_three, _d2)));
            // set state labels
            _m2 = _andnot(_three, _m2);
            _i2 = _or(_andnot(_three, _i2), _srli_si<1>(_three));
            _d2 = _or(_andnot(_three, _d2), _three);
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
        return (minscore + bitmask) >> trace_bits;
    }
};

} // namespace simd
} // namespace hmm
} // namespace octopus

#endif
