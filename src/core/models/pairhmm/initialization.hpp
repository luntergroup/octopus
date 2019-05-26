// Copyright (c) 2015-2019 Daniel Cooke and Gerton Lunter
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef initialization_hpp
#define initialization_hpp

#if __GNUC__ >= 6
    #pragma GCC diagnostic ignored "-Wignored-attributes"
#endif

namespace octopus { namespace hmm { namespace simd {


// Original rolling initialization using shift instructions
template<typename Policy>
class ShiftingRollingInit : Policy {
    using ScoreType  = typename Policy::ScoreType;
    using VectorType = typename Policy::VectorType;
    using Policy::vectorise_zero_set_last;
    using Policy::_left_shift_word;
    using Policy::_or;
    using Policy::_andnot;
    VectorType _initmask;
    VectorType _initmask2;
public:
    ShiftingRollingInit( ScoreType null_score_ ) {
        _initmask  = vectorise_zero_set_last(-1);
        _initmask2 = vectorise_zero_set_last(null_score_);
    }
    VectorType init( const VectorType& lhs ) const {
        return _or(_initmask2, _andnot(_initmask, lhs));
    }
    void update() {
        _initmask = _left_shift_word(_initmask);
        _initmask2 = _left_shift_word(_initmask2);
    }
};


// Faster rolling initialization using insert instructions
template<typename Policy>
class InsertRollingInit : Policy {
    using ScoreType  = typename Policy::ScoreType;
    using VectorType = typename Policy::VectorType;
    using Policy::_insert;
    int i;
    ScoreType value;
public:
    InsertRollingInit( ScoreType null_score_) : i {0}, value {null_score_} {}
    VectorType init( const VectorType& lhs ) const {
        return _insert(lhs, value, i);
    }
    void update() {
        ++i;
    }
};


} // simd
} // hmm
} // octopus

#endif
