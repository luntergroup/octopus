// Copyright (c) 2015-2021 Daniel Cooke and Gerton Lunter
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
#include <limits>
#include <vector>

#include <boost/align/aligned_allocator.hpp>

namespace octopus { namespace hmm { namespace simd {

template <typename InstructionSet,
          template <class> class InitializerType>
class PairHMM : private InstructionSet
{
public:
    using ScoreType = typename InstructionSet::ScoreType;

private:
    using VectorType  = typename InstructionSet::VectorType;
    using Initializer = InitializerType<InstructionSet>;
    using SmallVector = std::vector<VectorType, boost::alignment::aligned_allocator<VectorType>>;
    
    struct NullType {};
    
    // Methods
    using InstructionSet::vectorise;
    using InstructionSet::vectorise_zero_set_last;
    using InstructionSet::_extract;
    using InstructionSet::_insert_bottom;
    using InstructionSet::_insert_top;
    using InstructionSet::_add;
    using InstructionSet::_and;
    using InstructionSet::_andnot;
    using InstructionSet::_or;
    using InstructionSet::_cmpeq;
    using InstructionSet::_min;
    using InstructionSet::_max;
    using InstructionSet::_left_shift_word;
    using InstructionSet::_right_shift_word;
    
    // Constants
    constexpr static const char* name_ = InstructionSet::name;
    constexpr static int band_size_ {InstructionSet::band_size};
    constexpr static ScoreType infinity_tolerance_ {0x7FF};
    constexpr static ScoreType infinity_ {std::numeric_limits<ScoreType>::max() - infinity_tolerance_};
    constexpr static int trace_bits_ {2};
    constexpr static ScoreType n_score_ {2 << trace_bits_};
    
    constexpr static ScoreType max_quality_score_ {64}; // maximum reasonable phred base quality
    constexpr static ScoreType null_score_ {std::numeric_limits<ScoreType>::min()}; // baseline for score zero
    
    constexpr static int match_label_  {0};
    constexpr static int insert_label_ {1};
    constexpr static int delete_label_ {3};
    
    constexpr auto vectorise(NullType) const noexcept { return NullType {}; }
    
    template <int n, typename T>
    static auto vectorise_left_shift_bits(const T* values) noexcept
    {
        return _left_shift_bits<n>(vectorise(values));
    }
    template <int n>
    static auto vectorise_left_shift_bits(const std::int8_t value) noexcept
    {
        return vectorise(value << n);
    }
    template <int n>
    constexpr auto vectorise_left_shift_bits(NullType) const noexcept { return NullType {}; }
    template <int idx>
    static auto _right_shift_bits(const VectorType& vec) noexcept
    {
        return InstructionSet::template _right_shift_bits<idx>(vec);
    }
    template <int idx>
    static auto _left_shift_bits(const VectorType& vec) noexcept
    {
        return InstructionSet::template _left_shift_bits<idx>(vec);
    }
    
    VectorType maybe_right_shift_word(const VectorType& vec, const std::int8_t*) const noexcept
    {
        return _right_shift_word(vec);
    }
    const VectorType& maybe_right_shift_word(const VectorType& vec, ScoreType) const noexcept
    {
        return vec;
    }

    void update_gap_penalty(VectorType& current, const std::int8_t* source, const std::size_t gap_idx) const noexcept
    {
        current = _insert_top(_right_shift_word(current), source[gap_idx] << trace_bits_);
    }
    void update_gap_penalty(VectorType& current, const ScoreType source, const std::size_t gap_idx) const noexcept {}
    
    void
    update_snv_mask(VectorType& curr_mask,
                    VectorType& curr_caps,
                    const char* snv_mask,
                    const std::int8_t* caps,
                    const int pos,
                    const int truth_len) const noexcept
    {
        const bool pos_in_range {pos < truth_len};
        curr_mask = _insert_top(_right_shift_word(curr_mask), pos_in_range ? snv_mask[pos] : 'N');
        curr_caps = _insert_top(_right_shift_word(curr_caps), (pos_in_range ? caps[pos] : infinity_) << trace_bits_);
    }
    void update_snv_mask(NullType, NullType, NullType, NullType, int, int) const noexcept {}
    
    void
    update_match_state(VectorType& current,
                       const VectorType& _targetwin,
                       const VectorType& _truthwin,
                       const VectorType& _qualitieswin,
                       const VectorType& _truthnqual,
                       const VectorType& _snvmaskwin,
                       const VectorType& _snv_priorwin) const noexcept
    {
        auto _snvmask = _cmpeq(_targetwin, _snvmaskwin);
        current = _add(current, _min(_andnot(_cmpeq(_targetwin, _truthwin), _min(_qualitieswin, _or(_and(_snvmask, _snv_priorwin), _andnot(_snvmask, _qualitieswin)))), _truthnqual));
    }
    void
    update_match_state(VectorType& current,
                       const VectorType& _targetwin,
                       const VectorType& _truthwin,
                       const VectorType& _qualitieswin,
                       const VectorType& _truthnqual,
                       NullType, NullType) const noexcept
    {
        current = _add(current, _min(_andnot(_cmpeq(_targetwin, _truthwin), _qualitieswin), _truthnqual));
    }
    
    auto make_traceback_array(int target_len, int) const noexcept { return SmallVector(2 * (target_len + band_size_) + 1); }
    auto make_traceback_array(int target_len, NullType) const noexcept { return NullType {}; }
    
    void
    update_traceback(SmallVector& _backpointers,
                     const int index,
                     VectorType& match_state,
                     VectorType& insert_state,
                     VectorType& delete_state) const noexcept
    {
        const static VectorType _three = vectorise(3);
        const static VectorType _one = _right_shift_bits<1>(_three);   // could save one register
        _backpointers[index] = _or(_or(_and(_three, match_state), _left_shift_bits<2 * insert_label_>(_and(_three, insert_state))),
                               _left_shift_bits<2 * delete_label_>(_and(_three, delete_state)));
        // set state labels
        match_state = _andnot(_three, match_state);
        insert_state = _or(_andnot(_three, insert_state), _one);
        delete_state = _or(_andnot(_three, delete_state), _three);
    }
    void update_traceback(NullType, int, VectorType&, VectorType&, VectorType&) const noexcept {}
    
    void
    set_alignments(const char* truth,
                   const char* target,
                   const int truth_len,
                   const int target_len,
                   SmallVector& _backpointers,
                   const int minscoreidx,
                   int& first_pos,
                   char* align1,
                   char* align2) const noexcept
    {
        if (minscoreidx < 0) {
            // minscore was never updated so we must have overflowed badly
            first_pos = -1;
        }
        auto s = minscoreidx;    // point to the dummy match transition
        auto i      = s / 2 - target_len;
        auto y      = target_len;
        auto x      = s - y;
        auto alnidx = 0;
        auto ptr = reinterpret_cast<ScoreType*>(_backpointers.data() + s);
        if ((ptr + i) < reinterpret_cast<ScoreType*>(_backpointers.data())
            || (ptr + i) >= reinterpret_cast<ScoreType*>(_backpointers.data() + _backpointers.size())) {
            first_pos = -1;
            return;
        }
        auto state = (ptr[i] >> (2 * match_label_)) & 3;
        s -= 2;
        // this is 2*y (s even) or 2*y+1 (s odd)
        while (y > 0) {
            if (s < 0 || i < 0) {
                // This should never happen so must have overflowed
                first_pos = -1;
                return;
            }
            const auto new_state = (reinterpret_cast<ScoreType*>(_backpointers.data() + s)[i] >> (2 * state)) & 3;
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
    }
    void set_alignments(const char*, const char*, int, int, NullType, int, NullType, NullType, NullType) const noexcept {}
    
    template <typename OpenPenaltyArrayOrConstant,
              typename ExtendPenaltyArrayOrConstant,
              typename SnvMaskArrayOrNull,
              typename SnvBaseQualityCapArrayOrNull,
              typename PositionOrNull,
              typename CharArrayOrNull>
    int
    align_helper(const char* truth,
                 const char* target,
                 const std::int8_t* qualities,
                 const int truth_len,
                 const int target_len,
                 const OpenPenaltyArrayOrConstant gap_open,
                 const ExtendPenaltyArrayOrConstant gap_extend,
                 const SnvMaskArrayOrNull snv_mask,
                 const SnvBaseQualityCapArrayOrNull snv_prior,
                 const ScoreType nuc_prior,
                 PositionOrNull& first_pos,
                 CharArrayOrNull align1,
                 CharArrayOrNull align2) const noexcept
    {
        assert(target_len > 0 && truth_len > band_size_ && (truth_len == target_len + 2 * band_size_ - 1));
        const static VectorType _inf = vectorise(infinity_);
        const auto _nuc_prior = vectorise_left_shift_bits<trace_bits_>(nuc_prior);
        auto _truthwin     = vectorise(truth);
        auto _targetwin    = _inf;
        auto _qualitieswin = vectorise_left_shift_bits<trace_bits_>(max_quality_score_);
        auto _gap_open     = vectorise_left_shift_bits<trace_bits_>(gap_open);
        auto _gap_extend   = vectorise_left_shift_bits<trace_bits_>(gap_extend);
        auto _snvmaskwin   = vectorise(snv_mask);
        auto _snv_priorwin = vectorise_left_shift_bits<trace_bits_>(snv_prior);
        auto _truthnqual   = _add(_and(_cmpeq(_truthwin, vectorise('N')), vectorise(n_score_ - infinity_)), _inf);
        auto _backpointers = make_traceback_array(target_len, first_pos);
        auto _m1 = _inf, _i1 = _inf, _d1 = _inf, _m2 = _inf, _i2 = _inf, _d2 = _inf;
        Initializer rollinginit {null_score_};
        ScoreType minscore {infinity_}, cur_score;
        int minscoreidx {-1};
        for (int s = 0; s < 2 * (target_len + band_size_); s += 2) {
            // s even. truth is current; target needs updating
            _targetwin    = _left_shift_word(_targetwin);
            _qualitieswin = _left_shift_word(_qualitieswin);
            if (s / 2 < target_len) {
                _targetwin    = _insert_bottom(_targetwin, target[s / 2]);
                _qualitieswin = _insert_bottom(_qualitieswin, qualities[s / 2] << trace_bits_);
            } else {
                _targetwin    = _insert_bottom(_targetwin, '0');
                _qualitieswin = _insert_bottom(_qualitieswin, max_quality_score_ << trace_bits_);
            }
            _m1 = rollinginit.init(_m1);
            _m2 = rollinginit.init(_m2);
            _m1 = _min(_m1, _min(_i1, _d1));
            if (s / 2 >= target_len) {
                cur_score = _extract(_m1, s / 2 - target_len);
                if (cur_score < minscore) {
                    minscore = cur_score;
                    minscoreidx = s;
                }
            }
            update_match_state(_m1, _targetwin, _truthwin, _qualitieswin, _truthnqual, _snvmaskwin, _snv_priorwin);
            _d1 = _min(_add(_d2, maybe_right_shift_word(_gap_extend, gap_extend)), _add(_min(_m2, _i2), maybe_right_shift_word(_gap_open, gap_open))); // allow I->D
            _d1 = _insert_bottom(_left_shift_word(_d1), infinity_);
            _i1 = _add(_min(_add(_i2, _gap_extend), _add(_m2, _gap_open)), _nuc_prior);
            update_traceback(_backpointers, s, _m1, _i1, _d1);
            // S odd. Truth needs updating; target is current
            const auto pos = band_size_ + s / 2;
            const bool pos_in_range {pos < truth_len};
            const char base {pos_in_range ? truth[pos] : 'N'};
            _truthwin   = _insert_top(_right_shift_word(_truthwin), base);
            _truthnqual = _insert_top(_right_shift_word(_truthnqual), base == 'N' ? n_score_ : infinity_);
            const auto gap_idx = pos_in_range ? pos : truth_len - 1;
            update_gap_penalty(_gap_open, gap_open, gap_idx);
            update_gap_penalty(_gap_extend, gap_extend, gap_idx);
            update_snv_mask(_snvmaskwin, _snv_priorwin, snv_mask, snv_prior, pos, truth_len);
            rollinginit.update();
            _m2 = _min(_m2, _min(_i2, _d2));
            if (s / 2 >= target_len) {
                cur_score = _extract(_m2, s / 2 - target_len);
                if (cur_score < minscore) {
                    minscore = cur_score;
                    minscoreidx = s + 1;
                }
            }
            update_match_state(_m2, _targetwin, _truthwin, _qualitieswin, _truthnqual, _snvmaskwin, _snv_priorwin);
            _d2 = _min(_add(_d1, _gap_extend), _add(_min(_m1, _i1), _gap_open)); // allow I->D
            _i2 = _insert_top(_add(_min(_add(_right_shift_word(_i1), _gap_extend),
                                        _add(_right_shift_word(_m1), _gap_open)), _nuc_prior), infinity_);
            update_traceback(_backpointers, s + 1, _m2, _i2, _d2);
        }
        set_alignments(truth, target, truth_len, target_len, _backpointers, minscoreidx, first_pos, align1, align2);
        return (minscore - null_score_) >> trace_bits_;
    }
    
    auto get(const std::int8_t* values, int index) const noexcept { return values[index]; }
    auto get(std::int8_t value, int index) const noexcept { return value; }
    auto
    get_mismatch_quality(const char* target,
                         const std::int8_t* quals,
                         int x, int y,
                         const char* snv_mask,
                         const std::int8_t* caps) const noexcept
    {
        return (snv_mask[x] == target[y]) ? std::min(quals[y], caps[x]) : quals[y];
    }
    auto
    get_mismatch_quality(const char* target,
                         const std::int8_t* quals,
                         int x, int y,
                         NullType,
                         NullType) const noexcept
    {
        return quals[y];
    }
    
    template <typename CharArrayOrNull,
              typename OpenPenaltyArrayOrConstant,
              typename ExtendPenaltyArrayOrConstant,
              typename SnvMaskArrayOrNull,
              typename SnvBaseQualityCapArrayOrNull>
    int
    calculate_flank_score_helper(const int truth_len,
                                 const int lhs_flank_len,
                                 const int rhs_flank_len,
                                 const CharArrayOrNull target,
                                 const std::int8_t* quals,
                                 const OpenPenaltyArrayOrConstant gap_open,
                                 const ExtendPenaltyArrayOrConstant gap_extend,
                                 const SnvMaskArrayOrNull snv_mask,
                                 const SnvBaseQualityCapArrayOrNull snv_prior,
                                 const ScoreType nuc_prior,
                                 const int first_pos,
                                 const char* alignment1,
                                 const char* alignment2,
                                 int& target_mask_size) const noexcept
    {
        static constexpr char match {'M'}, insertion {'I'}, deletion {'D'};
        auto prev_state = match;
        int truth_idx {first_pos}, target_idx {0}, alignment_idx {0};
        int result {0};    // alignment score (within flank)
        const auto rhs_flank_begin = truth_len - rhs_flank_len;
        target_mask_size = 0;
        while (alignment1[alignment_idx]) {
            auto new_state = match;
            if (alignment1[alignment_idx] == gap_label) {
                new_state = insertion;
            } else if (alignment2[alignment_idx] == gap_label) { // can't be both '-'
                new_state = deletion;
            }
            const bool in_flank {truth_idx < lhs_flank_len || truth_idx >= rhs_flank_begin};
            switch (new_state) {
                case match:
                {
                    if (in_flank) {
                        if (alignment1[alignment_idx] != alignment2[alignment_idx]) {
                            if (alignment1[alignment_idx] != 'N') {
                                result += get_mismatch_quality(target, quals, truth_idx, target_idx, snv_mask, snv_prior);
                            } else {
                                result += n_score_ >> trace_bits_;
                            }
                        }
                        ++target_mask_size;
                    }
                    ++truth_idx;
                    ++target_idx;
                    break;
                }
                case insertion:
                {
                    if (in_flank) {
                        if (prev_state == insertion) {
                            result += get(gap_extend, truth_idx - 1) + nuc_prior;
                        } else {
                            // gap open score is charged for insertions just after the corresponding base, hence the -1
                            result += get(gap_open, truth_idx - 1) + nuc_prior;
                        }
                        ++target_mask_size;
                    }
                    ++target_idx;
                    break;
                }
                case deletion:
                {
                    if (in_flank) {
                        if (prev_state == deletion) {
                            result += get(gap_extend, truth_idx);
                        } else {
                            result += get(gap_open, truth_idx);
                        }
                    }
                    ++truth_idx;
                    break;
                }
            }
            ++alignment_idx;
            prev_state = new_state;
        }
        return result;
    }
    
public:
    constexpr static char gap_label = '-';
    
    constexpr static const char* name() noexcept { return name_; }
    constexpr static int band_size() noexcept { return band_size_; }
    
    template <typename OpenPenaltyArrayOrConstant,
              typename ExtendPenaltyArrayOrConstant>
    int
    align(const char* truth,
          const char* target,
          const std::int8_t* qualities,
          const int truth_len,
          const int target_len,
          const OpenPenaltyArrayOrConstant gap_open,
          const ExtendPenaltyArrayOrConstant gap_extend,
          short nuc_prior) const noexcept
    {
        constexpr static NullType null {};
        return align_helper(truth, target, qualities, truth_len, target_len, gap_open, gap_extend, null, null, nuc_prior, null, null, null);
    }
    
    template <typename OpenPenaltyArrayOrConstant,
              typename ExtendPenaltyArrayOrConstant>
    int
    align(const char* truth,
          const char* target,
          const std::int8_t* qualities,
          const int truth_len,
          const int target_len,
          const char* snv_mask,
          const std::int8_t* snv_prior,
          const OpenPenaltyArrayOrConstant gap_open,
          const ExtendPenaltyArrayOrConstant gap_extend,
          short nuc_prior) const noexcept
    {
        constexpr static NullType null {};
        return align_helper(truth, target, qualities, truth_len, target_len, gap_open, gap_extend, snv_mask, snv_prior, nuc_prior, null, null, null);
    }
    
    template <typename OpenPenaltyArrayOrConstant,
              typename ExtendPenaltyArrayOrConstant>
    int
    align(const char* truth,
          const char* target,
          const std::int8_t* qualities,
          const int truth_len,
          const int target_len,
          const OpenPenaltyArrayOrConstant gap_open,
          const ExtendPenaltyArrayOrConstant gap_extend,
          short nuc_prior,
          int& first_pos,
          char* align1,
          char* align2) const noexcept
    {
        constexpr static NullType null {};
        return align_helper(truth, target, qualities, truth_len, target_len, gap_open, gap_extend, null, null, nuc_prior, first_pos, align1, align2);
    }
    
    template <typename OpenPenaltyArrayOrConstant,
              typename ExtendPenaltyArrayOrConstant>
    int
    align(const char* truth,
          const char* target,
          const std::int8_t* qualities,
          const int truth_len,
          const int target_len,
          const char* snv_mask,
          const std::int8_t* snv_prior,
          const OpenPenaltyArrayOrConstant gap_open,
          const ExtendPenaltyArrayOrConstant gap_extend,
          short nuc_prior,
          int& first_pos,
          char* align1,
          char* align2) const noexcept
    {
        return align_helper(truth, target, qualities, truth_len, target_len, gap_open, gap_extend, snv_mask, snv_prior, nuc_prior, first_pos, align1, align2);
    }
    
    template <typename OpenPenaltyArrayOrConstant,
              typename ExtendPenaltyArrayOrConstant>
    int
    calculate_flank_score(const int truth_len,
                          const int lhs_flank_len,
                          const int rhs_flank_len,
                          const std::int8_t* quals,
                          const OpenPenaltyArrayOrConstant gap_open,
                          const ExtendPenaltyArrayOrConstant gap_extend,
                          const short nuc_prior,
                          const int first_pos,
                          const char* aln1,
                          const char* aln2,
                          int& target_mask_size) const noexcept
    {
        constexpr static NullType null {};
        return calculate_flank_score_helper(truth_len, lhs_flank_len, rhs_flank_len, null, quals, gap_open, gap_extend, null, null, nuc_prior, first_pos, aln1, aln2, target_mask_size);
    }
    
    template <typename OpenPenaltyArrayOrConstant,
              typename ExtendPenaltyArrayOrConstant>
    int
    calculate_flank_score(int truth_len,
                          int lhs_flank_len,
                          int rhs_flank_len,
                          const char* target,
                          const std::int8_t* quals,
                          const char* snv_mask,
                          const std::int8_t* snv_prior,
                          const OpenPenaltyArrayOrConstant gap_open,
                          const ExtendPenaltyArrayOrConstant gap_extend,
                          short nuc_prior,
                          int first_pos,
                          const char* aln1,
                          const char* aln2,
                          int& target_mask_size) const noexcept
    {
        return calculate_flank_score_helper(truth_len, lhs_flank_len, rhs_flank_len, target, quals, gap_open, gap_extend, snv_mask, snv_prior, nuc_prior, first_pos, aln1, aln2, target_mask_size);
    }
};

} // namespace simd
} // namespace hmm
} // namespace octopus

#endif
