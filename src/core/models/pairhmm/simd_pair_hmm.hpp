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
    
    using SmallVector = boost::container::small_vector<VectorType, 10000>;
    
    // Methods
    using InstructionSetPolicy::vectorise;
    using InstructionSetPolicy::vectorise_zero_set_last;
    using InstructionSetPolicy::_extract;
    using InstructionSetPolicy::_insert_bottom;
    using InstructionSetPolicy::_insert_top;
    using InstructionSetPolicy::_add;
    using InstructionSetPolicy::_and;
    using InstructionSetPolicy::_andnot;
    using InstructionSetPolicy::_or;
    using InstructionSetPolicy::_cmpeq;
    using InstructionSetPolicy::_min;
    using InstructionSetPolicy::_max;
    using InstructionSetPolicy::_left_shift_word;
    using InstructionSetPolicy::_right_shift_word;
    
    // Constants
    constexpr static int band_size_ = InstructionSetPolicy::band_size_;
    constexpr static ScoreType infinity_ = 0x7800;
    constexpr static int trace_bits_ = 2;
    constexpr static ScoreType n_score_ = 2 << trace_bits_;
    
    constexpr static int max_quality_score_ {64}; // maximum reasonable phred base quality
    constexpr static int null_score_ {-0x8000};   // baseline for score zero (== MIN_SHORT)
    
    constexpr static int match_label_  {0};
    constexpr static int insert_label_ {1};
    constexpr static int delete_label_ {3};
    
    const VectorType _inf = vectorise(infinity_);
    const VectorType _nscore_m_inf = vectorise(n_score_ - infinity_);
    const VectorType _n = vectorise('N');
    const VectorType _three = vectorise(3);
    const VectorType _one = _right_shift_bits<1>(_three);   // could save one register
    
    template <int n, typename T>
    auto vectorise_left_shift_bits(const T* values) const noexcept
    {
        return _left_shift_bits<n>(vectorise(values));
    }
    template <int n>
    auto vectorise_left_shift_bits(const std::int8_t value) const noexcept
    {
        return vectorise(value << n);
    }
    template <int idx>
    auto _right_shift_bits(const VectorType& vec) const noexcept {
        return InstructionSetPolicy::template _right_shift_bits<idx>(vec);
    }
    template <int idx>
    auto _left_shift_bits(const VectorType& vec) const noexcept {
        return InstructionSetPolicy::template _left_shift_bits<idx>(vec);
    }
    
    void update_gap_penalty(VectorType& current, const std::int8_t* source, const std::size_t gap_idx) const noexcept
    {
        current = _insert_top(_right_shift_word(current), source[gap_idx] << trace_bits_);
    }
    void update_gap_penalty(VectorType& current, const short source, const std::size_t gap_idx) const noexcept {}

public:
    constexpr static char gap_label = '-';
    
    constexpr int band_size() const noexcept { return band_size_; }
    
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
        assert(truth_len > band_size_ && (truth_len == target_len + 2 * band_size_ - 1));
        auto _m1 = _inf, _i1 = _inf, _d1 = _inf, _m2 = _inf, _i2 = _inf, _d2 = _inf;
        const auto _nuc_prior  = vectorise_left_shift_bits<trace_bits_>(nuc_prior);
        auto _initmask     = vectorise_zero_set_last(-1);
        auto _initmask2    = vectorise_zero_set_last(null_score_);
        auto _truthwin     = vectorise(truth);
        auto _targetwin    = _inf;
        auto _qualitieswin = vectorise_left_shift_bits<trace_bits_>(max_quality_score_);
        auto _gap_open     = vectorise_left_shift_bits<trace_bits_>(gap_open);
        auto _gap_extend   = vectorise_left_shift_bits<trace_bits_>(gap_extend);
        auto _truthnqual   = _add(_and(_cmpeq(_truthwin, _n), _nscore_m_inf), _inf);
        ScoreType minscore {infinity_};
        for (int s {0}; s <= 2 * (target_len + band_size_); s += 2) {
            // truth is current; target needs updating
            _targetwin    = _left_shift_word(_targetwin);
            _qualitieswin = _left_shift_word(_qualitieswin);
            if (s / 2 < target_len) {
                _targetwin    = _insert_bottom(_targetwin, target[s / 2]);
                _qualitieswin = _insert_bottom(_qualitieswin, qualities[s / 2] << trace_bits_);
            } else {
                _targetwin    = _insert_bottom(_targetwin, '0');
                _qualitieswin = _insert_bottom(_qualitieswin, max_quality_score_ << trace_bits_);
            }
            // S even
            _m1 = _or(_initmask2, _andnot(_initmask, _m1));
            _m2 = _or(_initmask2, _andnot(_initmask, _m2));
            _m1 = _min(_m1, _min(_i1, _d1));
            if (s / 2 >= target_len) {
                minscore = std::min(static_cast<decltype(minscore)>(_extract(_m1, std::max(0, s / 2 - target_len))), minscore);
            }
            _m1 = _add(_m1, _min(_andnot(_cmpeq(_targetwin, _truthwin), _qualitieswin), _truthnqual));
            _d1 = _min(_add(_d2, _gap_extend), _add(_min(_m2, _i2), _right_shift_word(_gap_open))); // allow I->D
            _d1 = _insert_bottom(_left_shift_word(_d1), infinity_);
            _i1 = _add(_min(_add(_i2, _gap_extend), _add(_m2, _gap_open)), _nuc_prior);
            // S odd
            // truth needs updating; target is current
            const auto pos = band_size_ + s / 2;
            const bool pos_in_range {pos < truth_len};
            const char base {pos_in_range ? truth[pos] : 'N'};
            _truthwin   = _insert_top(_right_shift_word(_truthwin), base);
            _truthnqual = _insert_top(_right_shift_word(_truthnqual), base == 'N' ? n_score_ : infinity_);
            const auto gap_idx = pos_in_range ? pos : truth_len - 1;
            update_gap_penalty(_gap_open, gap_open, gap_idx);
            update_gap_penalty(_gap_extend, gap_extend, gap_idx);
            _initmask  = _left_shift_word(_initmask);
            _initmask2 = _left_shift_word(_initmask2);
            _m2 = _min(_m2, _min(_i2, _d2));
            if (s / 2 >= target_len) {
                minscore = std::min(static_cast<decltype(minscore)>(_extract(_m2, s / 2 - target_len)), minscore);
            }
            _m2 = _add(_m2, _min(_andnot(_cmpeq(_targetwin, _truthwin), _qualitieswin), _truthnqual));
            _d2 = _min(_add(_d1, _gap_extend), _add(_min(_m1, _i1), _gap_open)); // allow I->D
            _i2 = _insert_top(_add(_min(_add(_right_shift_word(_i1), _gap_extend),
                                        _add(_right_shift_word(_m1), _gap_open)), _nuc_prior), infinity_);
        }
        return (minscore - null_score_) >> trace_bits_;
    }
    
    template <typename OpenPenalty,
              typename ExtendPenalty>
    int align(const char* truth,
              const char* target,
              const std::int8_t* qualities,
              const int truth_len,
              const int target_len,
              const char* snv_mask,
              const std::int8_t* snv_prior,
              const OpenPenalty gap_open,
              const ExtendPenalty gap_extend,
              short nuc_prior) const noexcept
    {
        assert(truth_len > band_size_ && (truth_len == target_len + 2 * band_size_ - 1));
        auto _m1 = _inf, _i1 = _inf, _d1 = _inf, _m2 = _inf, _i2 = _inf, _d2 = _inf;
        const auto _nuc_prior = vectorise_left_shift_bits<trace_bits_>(nuc_prior);
        auto _initmask     = vectorise_zero_set_last(-1);
        auto _initmask2    = vectorise_zero_set_last(null_score_);
        auto _truthwin     = vectorise(truth);
        auto _targetwin    = _inf;
        auto _qualitieswin = vectorise_left_shift_bits<trace_bits_>(max_quality_score_);
        auto _gap_open     = vectorise_left_shift_bits<trace_bits_>(gap_open);
        auto _gap_extend   = vectorise_left_shift_bits<trace_bits_>(gap_extend);
        auto _snvmaskwin   = vectorise(snv_mask);
        auto _snv_priorwin = vectorise_left_shift_bits<trace_bits_>(snv_prior);
        auto _truthnqual   = _add(_and(_cmpeq(_truthwin, _n), _nscore_m_inf), _inf);
        VectorType _snvmask;
        ScoreType minscore {infinity_};
        for (int s {0}; s <= 2 * (target_len + band_size_); s += 2) {
            // truth is current; target needs updating
            _targetwin    = _left_shift_word(_targetwin);
            _qualitieswin = _left_shift_word(_qualitieswin);
            if (s / 2 < target_len) {
                _targetwin    = _insert_bottom(_targetwin, target[s / 2]);
                _qualitieswin = _insert_bottom(_qualitieswin, qualities[s / 2] << trace_bits_);
            } else {
                _targetwin    = _insert_bottom(_targetwin, '0');
                _qualitieswin = _insert_bottom(_qualitieswin, max_quality_score_ << trace_bits_);
            }
            // S even
            _m1 = _or(_initmask2, _andnot(_initmask, _m1));
            _m2 = _or(_initmask2, _andnot(_initmask, _m2));
            _m1 = _min(_m1, _min(_i1, _d1));
            if (s / 2 >= target_len) {
                minscore = std::min(static_cast<decltype(minscore)>(_extract(_m1, std::max(0, s / 2 - target_len))), minscore);
            }
            _snvmask = _cmpeq(_targetwin, _snvmaskwin);
            _m1 = _add(_m1, _min(_andnot(_cmpeq(_targetwin, _truthwin), _min(_qualitieswin, _or(_and(_snvmask, _snv_priorwin), _andnot(_snvmask, _qualitieswin)))), _truthnqual));
            _d1 = _min(_add(_d2, _gap_extend), _add(_min(_m2, _i2), _right_shift_word(_gap_open))); // allow I->D
            _d1 = _insert_bottom(_left_shift_word(_d1), infinity_);
            _i1 = _add(_min(_add(_i2, _gap_extend), _add(_m2, _gap_open)), _nuc_prior);
            // S odd
            // truth needs updating; target is current
            const auto pos = band_size_ + s / 2;
            const bool pos_in_range {pos < truth_len};
            const char base {pos_in_range ? truth[pos] : 'N'};
            _truthwin   = _insert_top(_right_shift_word(_truthwin), base);
            _truthnqual = _insert_top(_right_shift_word(_truthnqual), base == 'N' ? n_score_ : infinity_);
            const auto gap_idx = pos_in_range ? pos : truth_len - 1;
            update_gap_penalty(_gap_open, gap_open, gap_idx);
            update_gap_penalty(_gap_extend, gap_extend, gap_idx);
            _snvmaskwin   = _insert_top(_right_shift_word(_snvmaskwin), pos_in_range ? snv_mask[pos] : 'N');
            _snv_priorwin = _insert_top(_right_shift_word(_snv_priorwin), (pos_in_range ? snv_prior[pos] : infinity_) << trace_bits_);
            _initmask  = _left_shift_word(_initmask);
            _initmask2 = _left_shift_word(_initmask2);
            _m2 = _min(_m2, _min(_i2, _d2));
            if (s / 2 >= target_len) {
                minscore = std::min(static_cast<decltype(minscore)>(_extract(_m2, s / 2 - target_len)), minscore);
            }
            _snvmask = _cmpeq(_targetwin, _snvmaskwin);
            _m2 = _add(_m2, _min(_andnot(_cmpeq(_targetwin, _truthwin), _min(_qualitieswin, _or(_and(_snvmask, _snv_priorwin), _andnot(_snvmask, _qualitieswin)))), _truthnqual));
            _d2 = _min(_add(_d1, _gap_extend), _add(_min(_m1, _i1), _gap_open)); // allow I->D
            _i2 = _insert_top(_add(_min(_add(_right_shift_word(_i1), _gap_extend),
                                        _add(_right_shift_word(_m1), _gap_open)), _nuc_prior), infinity_);
        }
        return (minscore - null_score_) >> trace_bits_;
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
          char* align2) const noexcept
    {
        assert(truth_len > band_size_ && (truth_len == target_len + 2 * band_size_ - 1));
        auto _m1 = _inf, _i1 = _inf, _d1 = _inf, _m2 = _inf, _i2 = _inf, _d2 = _inf;
        const auto _nuc_prior = vectorise_left_shift_bits<trace_bits_>(nuc_prior);
        auto _initmask     = vectorise_zero_set_last(-1);
        auto _initmask2    = vectorise_zero_set_last(null_score_);
        auto _truthwin     = vectorise(truth);
        auto _targetwin    = _inf;
        auto _qualitieswin = vectorise_left_shift_bits<trace_bits_>(max_quality_score_);
        auto _gap_open     = vectorise_left_shift_bits<trace_bits_>(gap_open);
        auto _gap_extend   = vectorise_left_shift_bits<trace_bits_>(gap_extend);
        auto _truthnqual   = _add(_and(_cmpeq(_truthwin, _n), _nscore_m_inf), _inf);
        SmallVector _backpointers(2 * (truth_len + band_size_));
        ScoreType minscore {infinity_}, cur_score;
        int s, minscoreidx {-1};
        for (s = 0; s <= 2 * (target_len + band_size_); s += 2) {
            // truth is current; target needs updating
            _targetwin    = _left_shift_word(_targetwin);
            _qualitieswin = _left_shift_word(_qualitieswin);
            if (s / 2 < target_len) {
                _targetwin    = _insert_bottom(_targetwin, target[s / 2]);
                _qualitieswin = _insert_bottom(_qualitieswin, qualities[s / 2] << trace_bits_);
            } else {
                _targetwin    = _insert_bottom(_targetwin, '0');
                _qualitieswin = _insert_bottom(_qualitieswin, max_quality_score_ << trace_bits_);
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
            _d1 = _min(_add(_d2, _gap_extend), _add(_min(_m2, _i2), _right_shift_word(_gap_open))); // allow I->D
            _d1 = _insert_bottom(_left_shift_word(_d1), infinity_);
            _i1 = _add(_min(_add(_i2, _gap_extend), _add(_m2, _gap_open)), _nuc_prior);
            _backpointers[s] = _or(_or(_and(_three, _m1), _left_shift_bits<2 * insert_label_>(_and(_three, _i1))),
                                   _left_shift_bits<2 * delete_label_>(_and(_three, _d1)));
            // set state labels
            _m1 = _andnot(_three, _m1);
            _i1 = _or(_andnot(_three, _i1), _one);
            _d1 = _or(_andnot(_three, _d1), _three);
            // S odd
            // truth needs updating; target is current
            const auto pos = band_size_ + s / 2;
            const bool pos_in_range {pos < truth_len};
            const char base {pos_in_range ? truth[pos] : 'N'};
            _truthwin   = _insert_top(_right_shift_word(_truthwin), base);
            _truthnqual = _insert_top(_right_shift_word(_truthnqual), base == 'N' ? n_score_ : infinity_);
            const auto gap_idx = pos_in_range ? pos : truth_len - 1;
            update_gap_penalty(_gap_open, gap_open, gap_idx);
            update_gap_penalty(_gap_extend, gap_extend, gap_idx);
            _initmask  = _left_shift_word(_initmask);
            _initmask2 = _left_shift_word(_initmask2);
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
            _i2 = _insert_top(_add(_min(_add(_right_shift_word(_i1), _gap_extend),
                                        _add(_right_shift_word(_m1), _gap_open)), _nuc_prior), infinity_);
            _backpointers[s + 1] = _or(_or(_and(_three, _m2), _left_shift_bits<2 * insert_label_>(_and(_three, _i2))),
                                       _left_shift_bits<2 * delete_label_>(_and(_three, _d2)));
            // set state labels
            _m2 = _andnot(_three, _m2);
            _i2 = _or(_andnot(_three, _i2), _one);
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
        auto state = (ptr[i] >> (2 * match_label_)) & 3;
        s -= 2;
        // this is 2*y (s even) or 2*y+1 (s odd)
        while (y > 0) {
            if (s < 0 || i < 0) {
                // This should never happen so must have overflowed
                first_pos = -1;
                return -1;
            }
            const auto new_state = (reinterpret_cast<short*>(_backpointers.data() + s)[i] >> (2 * state)) & 3;
            if (state == match_label_) {
                s -= 2;
                align1[alnidx] = truth[--x];
                align2[alnidx] = target[--y];
            } else if (state == insert_label_) {
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
        return (minscore - null_score_) >> trace_bits_;
    }
    
    template <typename OpenPenalty,
              typename ExtendPenalty>
    int
    align(const char* truth,
          const char* target,
          const std::int8_t* qualities,
          const int truth_len,
          const int target_len,
          const char* snv_mask,
          const std::int8_t* snv_prior,
          const OpenPenalty gap_open,
          const ExtendPenalty gap_extend,
          short nuc_prior,
          int& first_pos,
          char* align1,
          char* align2) const noexcept
    {
        assert(truth_len > band_size_ && (truth_len == target_len + 2 * band_size_ - 1));
        auto _m1 = _inf, _i1 = _inf, _d1 = _inf, _m2 = _inf, _i2 = _inf, _d2 = _inf;
        const auto _nuc_prior = vectorise_left_shift_bits<trace_bits_>(nuc_prior);
        auto _initmask     = vectorise_zero_set_last(-1);
        auto _initmask2    = vectorise_zero_set_last(null_score_);
        auto _truthwin     = vectorise(truth);
        auto _targetwin    = _inf;
        auto _qualitieswin = vectorise_left_shift_bits<trace_bits_>(max_quality_score_);
        auto _gap_open     = vectorise_left_shift_bits<trace_bits_>(gap_open);
        auto _gap_extend   = vectorise_left_shift_bits<trace_bits_>(gap_extend);
        auto _snvmaskwin   = vectorise(snv_mask);
        auto _snv_priorwin = vectorise_left_shift_bits<trace_bits_>(snv_prior);
        auto _truthnqual   = _add(_and(_cmpeq(_truthwin, _n), _nscore_m_inf), _inf);
        VectorType _snvmask;
        SmallVector _backpointers(2 * (truth_len + band_size_));
        ScoreType minscore {infinity_}, cur_score;
        int s, minscoreidx {-1};
        for (s = 0; s <= 2 * (target_len + band_size_); s += 2) {
            // truth is current; target needs updating
            _targetwin    = _left_shift_word(_targetwin);
            _qualitieswin = _left_shift_word(_qualitieswin);
            if (s / 2 < target_len) {
                _targetwin    = _insert_bottom(_targetwin, target[s / 2]);
                _qualitieswin = _insert_bottom(_qualitieswin, qualities[s / 2] << trace_bits_);
            } else {
                _targetwin    = _insert_bottom(_targetwin, '0');
                _qualitieswin = _insert_bottom(_qualitieswin, max_quality_score_ << trace_bits_);
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
            _snvmask = _cmpeq(_targetwin, _snvmaskwin);
            _m1 = _add(_m1, _min(_andnot(_cmpeq(_targetwin, _truthwin), _min(_qualitieswin, _or(_and(_snvmask, _snv_priorwin), _andnot(_snvmask, _qualitieswin)))), _truthnqual));
            _d1 = _min(_add(_d2, _gap_extend), _add(_min(_m2, _i2), _right_shift_word(_gap_open))); // allow I->D
            _d1 = _insert_bottom(_left_shift_word(_d1), infinity_);
            _i1 = _add(_min(_add(_i2, _gap_extend), _add(_m2, _gap_open)), _nuc_prior);
            _backpointers[s] = _or(_or(_and(_three, _m1), _left_shift_bits<2 * insert_label_>(_and(_three, _i1))),
                                   _left_shift_bits<2 * delete_label_>(_and(_three, _d1)));
            // set state labels
            _m1 = _andnot(_three, _m1);
            _i1 = _or(_andnot(_three, _i1), _one);
            _d1 = _or(_andnot(_three, _d1), _three);
            // S odd
            // truth needs updating; target is current
            const auto pos = band_size_ + s / 2;
            const bool pos_in_range {pos < truth_len};
            const char base {pos_in_range ? truth[pos] : 'N'};
            _truthwin   = _insert_top(_right_shift_word(_truthwin), base);
            _truthnqual = _insert_top(_right_shift_word(_truthnqual), base == 'N' ? n_score_ : infinity_);
            const auto gap_idx = pos_in_range ? pos : truth_len - 1;
            update_gap_penalty(_gap_open, gap_open, gap_idx);
            update_gap_penalty(_gap_extend, gap_extend, gap_idx);
            _snvmaskwin   = _insert_top(_right_shift_word(_snvmaskwin), pos_in_range ? snv_mask[pos] : 'N');
            _snv_priorwin = _insert_top(_right_shift_word(_snv_priorwin), (pos_in_range ? snv_prior[pos] : infinity_) << trace_bits_);
            _initmask  = _left_shift_word(_initmask);
            _initmask2 = _left_shift_word(_initmask2);
            _m2 = _min(_m2, _min(_i2, _d2));
            if (s / 2 >= target_len) {
                cur_score = _extract(_m2, s / 2 - target_len);
                if (cur_score < minscore) {
                    minscore = cur_score;
                    minscoreidx = s + 1;
                }
            }
            _snvmask = _cmpeq(_targetwin, _snvmaskwin);
            _m2 = _add(_m2, _min(_andnot(_cmpeq(_targetwin, _truthwin), _min(_qualitieswin, _or(_and(_snvmask, _snv_priorwin), _andnot(_snvmask, _qualitieswin)))), _truthnqual));
            _d2 = _min(_add(_d1, _gap_extend), _add(_min(_m1, _i1), _gap_open)); // allow I->D
            _i2 = _insert_top(_add(_min(_add(_right_shift_word(_i1), _gap_extend),
                                        _add(_right_shift_word(_m1), _gap_open)), _nuc_prior), infinity_);
            _backpointers[s + 1] = _or(_or(_and(_three, _m2), _left_shift_bits<2 * insert_label_>(_and(_three, _i2))),
                                       _left_shift_bits<2 * delete_label_>(_and(_three, _d2)));
            // set state labels
            _m2 = _andnot(_three, _m2);
            _i2 = _or(_andnot(_three, _i2), _one);
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
        auto state = (ptr[i] >> (2 * match_label_)) & 3;
        s -= 2;
        // this is 2*y (s even) or 2*y+1 (s odd)
        while (y > 0) {
            if (s < 0 || i < 0) {
                // This should never happen so must have overflowed
                first_pos = -1;
                return -1;
            }
            const auto new_state = (reinterpret_cast<short*>(_backpointers.data() + s)[i] >> (2 * state)) & 3;
            if (state == match_label_) {
                s -= 2;
                align1[alnidx] = truth[--x];
                align2[alnidx] = target[--y];
            } else if (state == insert_label_) {
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
        return (minscore - null_score_) >> trace_bits_;
    }
    
    int
    calculate_flank_score(int truth_len,
                          int lhs_flank_len,
                          int rhs_flank_len,
                          const char* target,
                          const std::int8_t* quals,
                          const char* snv_mask,
                          const std::int8_t* snv_prior,
                          const std::int8_t* gap_open,
                          const std::int8_t* gap_extend,
                          short nuc_prior,
                          int first_pos,
                          const char* aln1,
                          const char* aln2,
                          int& target_mask_size) const noexcept
    {
        static constexpr char match {'M'}, insertion {'I'}, deletion {'D'};
        auto prev_state = match;
        int x {first_pos}; // index into truth
        int y {0};         // index into target
        int i {0};         // index into alignment
        int result {0};    // alignment score (within flank)
        const auto rhs_flank_begin = truth_len - rhs_flank_len;
        target_mask_size = 0;
        while (aln1[i]) {
            auto new_state = match;
            if (aln1[i] == gap_label) {
                new_state = insertion;
            } else if (aln2[i] == gap_label) { // can't be both '-'
                new_state = deletion;
            }
            switch (new_state) {
                case match:
                {
                    if (x < lhs_flank_len || x >= rhs_flank_begin) {
                        if (aln1[i] != aln2[i]) {
                            if (aln1[i] != 'N') {
                                result += (snv_mask[x] == target[y]) ? std::min(quals[y], snv_prior[x]) : quals[y];
                            } else {
                                result += n_score_ >> trace_bits_;
                            }
                        }
                        ++target_mask_size;
                    }
                    ++x;
                    ++y;
                    break;
                }
                case insertion:
                {
                    if (x < lhs_flank_len || x >= rhs_flank_begin) {
                        if (prev_state == insertion) {
                            result += gap_extend[x - 1] + nuc_prior;
                        } else {
                            // gap open score is charged for insertions just after the corresponding base, hence the -1
                            result += gap_open[x - 1] + nuc_prior;
                        }
                        ++target_mask_size;
                    }
                    ++y;
                    break;
                }
                case deletion:
                {
                    if (x < lhs_flank_len || x >= rhs_flank_begin) {
                        if (prev_state == deletion) {
                            result += gap_extend[x];
                        } else {
                            result += gap_open[x];
                        }
                    }
                    ++x;
                    break;
                }
            }
            ++i;
            prev_state = new_state;
        }
        return result;
    }
    
    int
    calculate_flank_score(const int truth_len,
                          const int lhs_flank_len,
                          const int rhs_flank_len,
                          const std::int8_t* quals,
                          const std::int8_t* gap_open,
                          const short gap_extend,
                          const short nuc_prior,
                          const int first_pos,
                          const char* aln1,
                          const char* aln2,
                          int& target_mask_size) const noexcept
    {
        static constexpr char match {'M'}, insertion {'I'}, deletion {'D'};
        auto prev_state = match;
        int x {first_pos}; // index into haplotype
        int y {0};         // index into read
        int i {0};         // index into alignment
        int result {0};    // alignment score (within flank)
        target_mask_size = 0;
        while (aln1[i]) {
            auto new_state = match;
            if (aln1[i] == gap_label) {
                new_state = insertion;
            } else if (aln2[i] == gap_label) { // can't be both '-'
                new_state = deletion;
            }
            switch (new_state) {
                case match:
                {
                    if (x < lhs_flank_len || x >= (truth_len - rhs_flank_len)) {
                        if (aln1[i] != aln2[i]) {
                            if (aln1[i] != 'N') {
                                result += quals[y];
                            } else {
                                result += n_score_ >> trace_bits_;
                            }
                        }
                        ++target_mask_size;
                    }
                    ++x;
                    ++y;
                    break;
                }
                case insertion:
                {
                    if (x < lhs_flank_len || x >= (truth_len - rhs_flank_len)) {
                        if (prev_state == insertion) {
                            result += gap_extend + nuc_prior;
                        } else {
                            // gap open score is charged for insertions just after the corresponding base, hence the -1
                            result += gap_open[x - 1] + nuc_prior;
                        }
                        ++target_mask_size;
                    }
                    ++y;
                    break;
                }
                case deletion:
                {
                    if (x < lhs_flank_len || x >= (truth_len - rhs_flank_len)) {
                        if (prev_state == deletion) {
                            result += gap_extend;
                        } else {
                            result += gap_open[x];
                        }
                    }
                    ++x;
                    break;
                }
            }
            
            ++i;
            prev_state = new_state;
        }
        return result;
    }
    
    int
    calculate_flank_score(const int truth_len,
                          const int lhs_flank_len,
                          const int rhs_flank_len,
                          const char* target,
                          const std::int8_t* quals,
                          const char* snv_mask,
                          const std::int8_t* snv_prior,
                          const std::int8_t* gap_open,
                          const short gap_extend,
                          const short nuc_prior,
                          const int first_pos,
                          const char* aln1,
                          const char* aln2,
                          int& target_mask_size) const noexcept
    {
        static constexpr char match {'M'}, insertion {'I'}, deletion {'D'};
        auto prev_state = match;
        int x {first_pos}; // index into truth
        int y {0};         // index into target
        int i {0};         // index into alignment
        int result {0};    // alignment score (within flank)
        const auto rhs_flank_begin = truth_len - rhs_flank_len;
        target_mask_size = 0;
        while (aln1[i]) {
            auto new_state = match;
            if (aln1[i] == gap_label) {
                new_state = insertion;
            } else if (aln2[i] == gap_label) { // can't be both '-'
                new_state = deletion;
            }
            switch (new_state) {
                case match:
                {
                    if (x < lhs_flank_len || x >= rhs_flank_begin) {
                        if (aln1[i] != aln2[i]) {
                            if (aln1[i] != 'N') {
                                result += (snv_mask[x] == target[y]) ? std::min(quals[y], snv_prior[x]) : quals[y];
                            } else {
                                result += n_score_ >> trace_bits_;
                            }
                        }
                        ++target_mask_size;
                    }
                    ++x;
                    ++y;
                    break;
                }
                case insertion:
                {
                    if (x < lhs_flank_len || x >= rhs_flank_begin) {
                        if (prev_state == insertion) {
                            result += gap_extend + nuc_prior;
                        } else {
                            // gap open score is charged for insertions just after the corresponding base, hence the -1
                            result += gap_open[x - 1] + nuc_prior;
                        }
                        ++target_mask_size;
                    }
                    ++y;
                    break;
                }
                case deletion:
                {
                    if (x < lhs_flank_len || x >= rhs_flank_begin) {
                        if (prev_state == deletion) {
                            result += gap_extend;
                        } else {
                            result += gap_open[x];
                        }
                    }
                    ++x;
                    break;
                }
            }
            ++i;
            prev_state = new_state;
        }
        return result;
    }
};

} // namespace simd
} // namespace hmm
} // namespace octopus

#endif
