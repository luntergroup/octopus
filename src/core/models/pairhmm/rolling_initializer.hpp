// Copyright (c) 2015-2021 Daniel Cooke and Gerton Lunter
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef rolling_initializer_hpp
#define rolling_initializer_hpp

#if __GNUC__ >= 6
    #pragma GCC diagnostic ignored "-Wignored-attributes"
#endif

namespace octopus { namespace hmm { namespace simd {

// Original rolling initialization using shift instructions
template <typename PairHMMInstructionSet>
class ShiftingRollingInitializer : private PairHMMInstructionSet
{
    using ScoreType  = typename PairHMMInstructionSet::ScoreType;
    using VectorType = typename PairHMMInstructionSet::VectorType;
    using PairHMMInstructionSet::vectorise_zero_set_last;
    using PairHMMInstructionSet::_left_shift_word;
    using PairHMMInstructionSet::_or;
    using PairHMMInstructionSet::_andnot;
    VectorType _initmask_;
    VectorType _initmask2_;
public:
    ShiftingRollingInitializer(const ScoreType null_score) noexcept
    : _initmask_ {vectorise_zero_set_last(-1)}
    , _initmask2_ {vectorise_zero_set_last(null_score)}
    {}
    VectorType init(const VectorType& a) const noexcept { return _or(_initmask2_, _andnot(_initmask_, a)); }
    void update() noexcept
    {
        _initmask_  = _left_shift_word(_initmask_);
        _initmask2_ = _left_shift_word(_initmask2_);
    }
};

// Faster rolling initialization using insert instructions
template <typename PairHMMInstructionSet>
class InsertRollingInitializer : private PairHMMInstructionSet
{
    using ScoreType  = typename PairHMMInstructionSet::ScoreType;
    using VectorType = typename PairHMMInstructionSet::VectorType;
    using PairHMMInstructionSet::_insert;
    int i_;
    ScoreType value_;
public:
    constexpr InsertRollingInitializer(ScoreType null_score) noexcept : i_ {0}, value_ {null_score} {}
    VectorType init(const VectorType& a) const noexcept { return _insert(a, value_, i_); }
    void update() noexcept { ++i_; }
};

} // simd
} // hmm
} // octopus

#endif
