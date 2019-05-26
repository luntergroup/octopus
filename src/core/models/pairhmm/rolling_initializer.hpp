// Copyright (c) 2015-2019 Daniel Cooke and Gerton Lunter
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef rolling_initializer_hpp
#define rolling_initializer_hpp

#if __GNUC__ >= 6
    #pragma GCC diagnostic ignored "-Wignored-attributes"
#endif

namespace octopus { namespace hmm { namespace simd {

// Original rolling initialization using shift instructions
template <typename PairHMMInstructionSet>
class ShiftingRollingInitializer : PairHMMInstructionSet
{
    using ScoreType  = typename PairHMMInstructionSet::ScoreType;
    using VectorType = typename PairHMMInstructionSet::VectorType;
    using PairHMMInstructionSet::vectorise_zero_set_last;
    using PairHMMInstructionSet::_left_shift_word;
    using PairHMMInstructionSet::_or;
    using PairHMMInstructionSet::_andnot;
    VectorType _initmask;
    VectorType _initmask2;
public:
    ShiftingRollingInitializer(ScoreType null_score_) noexcept
    : _initmask {vectorise_zero_set_last(-1)}
    , _initmask2 {vectorise_zero_set_last(null_score_)}
    {}
    VectorType init(const VectorType& a) const noexcept { return _or(_initmask2, _andnot(_initmask, a)); }
    void update() noexcept
    {
        _initmask  = _left_shift_word(_initmask);
        _initmask2 = _left_shift_word(_initmask2);
    }
};

// Faster rolling initialization using insert instructions
template <typename PairHMMInstructionSet>
class InsertRollingInitializer : PairHMMInstructionSet
{
    using ScoreType  = typename PairHMMInstructionSet::ScoreType;
    using VectorType = typename PairHMMInstructionSet::VectorType;
    using PairHMMInstructionSet::_insert;
    int i;
    ScoreType value;
public:
    InsertRollingInitializer(ScoreType null_score_) noexcept : i {0}, value {null_score_} {}
    VectorType init(const VectorType& a) const noexcept { return _insert(a, value, i); }
    void update() noexcept { ++i; }
};

} // simd
} // hmm
} // octopus

#endif
