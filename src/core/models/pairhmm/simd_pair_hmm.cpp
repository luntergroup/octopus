// Copyright (c) 2017 Daniel Cooke and Gerton Lunter
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#if __GNUC__ >= 6
    #pragma GCC diagnostic ignored "-Wignored-attributes"
#endif

#include "simd_pair_hmm.hpp"

#include <vector>
#include <algorithm>
#include <emmintrin.h>
#include <cassert>

#include <boost/container/small_vector.hpp>

//#include <iostream> // DEBUG
//#include <iterator> // DEBUG
//
//struct Seq { __m128i val; };
//struct Qual { __m128i val; };
//struct Mask { __m128i val; };
//
//std::ostream& operator<<(std::ostream& os, const Seq s)
//{
//    const auto val = (std::uint16_t*) &s.val;
//    std::copy(val, val + 8, std::ostreambuf_iterator<char>(os));
//    return os;
//}
//
//std::ostream& operator<<(std::ostream& os, const Qual q)
//{
//    const auto val = (std::uint16_t*) &q.val;
//    std::transform(val, val + 8, std::ostream_iterator<unsigned>(os, " "),
//                   [] (const auto x) { return x >> 2; });
//    return os;
//}
//
//std::ostream& operator<<(std::ostream& os, const Mask m)
//{
//    const auto val = (std::uint16_t*) &m.val;
//    std::copy(val, val + 8, std::ostream_iterator<bool>(os));
//    return os;
//}

namespace octopus { namespace hmm { namespace simd {

constexpr std::size_t staticBackpointerCapacity {10000};

template <typename T>
using SmallVector = boost::container::small_vector<T, staticBackpointerCapacity>;

constexpr short nScore {2 << 2};
constexpr int bandSize {8};
constexpr short inf {0x7800};
constexpr char gap {'-'};

auto extract_epi16(const __m128i a, const int imm) noexcept
{
    switch (imm) {
        case 0:  return _mm_extract_epi16(a, 0);
        case 1:  return _mm_extract_epi16(a, 1);
        case 2:  return _mm_extract_epi16(a, 2);
        case 3:  return _mm_extract_epi16(a, 3);
        case 4:  return _mm_extract_epi16(a, 4);
        case 5:  return _mm_extract_epi16(a, 5);
        case 6:  return _mm_extract_epi16(a, 6);
        default: return _mm_extract_epi16(a, 7);
    }
}

int align(const char* truth, const char* target, const std::int8_t* qualities,
          const int truth_len, const int target_len,
          const std::int8_t* gap_open, short gap_extend, short nuc_prior) noexcept
{
    // target is the read; the shorter of the sequences
    // no checks for overflow are done
    //
    // the bottom-left and top-right corners of the DP table are just
    // included at the extreme ends of the diagonal, which measures
    // n=8 entries diagonally across.  This fixes the length of the
    // longer (horizontal) sequence to 15 (2*8-1) more than the shorter
    //
    // the << 2's are because the lower two bits are reserved for back tracing
    
    assert(truth_len > bandSize && (truth_len == target_len + 2 * bandSize - 1));
    
    gap_extend <<= 2;
    nuc_prior <<= 2;
    
    using SimdInt = __m128i;
    
    SimdInt _m1 {_mm_set1_epi16(inf)};
    auto _i1 = _m1;
    auto _d1 = _m1;
    auto _m2 = _m1;
    auto _i2 = _m1;
    auto _d2 = _m1;
    
    SimdInt _gap_extend {_mm_set1_epi16(gap_extend)};
    SimdInt _nuc_prior  {_mm_set1_epi16(nuc_prior)};
    
    SimdInt _initmask   {_mm_set_epi16(0,0,0,0,0,0,0,-1)};
    SimdInt _initmask2  {_mm_set_epi16(0,0,0,0,0,0,0,-0x8000)};
    
    // truth is initialized with the n-long prefix, in forward direction
    // target is initialized as empty; reverse direction
    SimdInt _truthwin  {_mm_set_epi16(truth[7], truth[6], truth[5], truth[4],
                                      truth[3], truth[2], truth[1], truth[0])};
    SimdInt _targetwin  {_m1};
    SimdInt _qualitieswin {_mm_set1_epi16(64 << 2)};
    
    // if N, make nScore; if != N, make inf
    SimdInt _truthnqual {_mm_add_epi16(_mm_and_si128(_mm_cmpeq_epi16(_truthwin, _mm_set1_epi16('N')),
                                                     _mm_set1_epi16(nScore - inf)),
                                       _mm_set1_epi16(inf))};
    
    SimdInt _gap_open {_mm_set_epi16(gap_open[7] << 2,gap_open[6] << 2,gap_open[5] << 2,gap_open[4] << 2,
                                     gap_open[3] << 2,gap_open[2] << 2,gap_open[1] << 2,gap_open[0] << 2)};
    
    short minscore {inf};
    
    for (int s {0}; s <= 2 * (target_len + bandSize); s += 2) {
        // truth is current; target needs updating
        _targetwin    = _mm_slli_si128(_targetwin, 2);
        _qualitieswin = _mm_slli_si128(_qualitieswin, 2);
        
        if (s / 2 < target_len) {
            _targetwin    = _mm_insert_epi16(_targetwin, target[s / 2], 0);
            _qualitieswin = _mm_insert_epi16(_qualitieswin, qualities[s / 2] << 2, 0);
        } else {
            _targetwin    = _mm_insert_epi16(_targetwin, '0', 0);
            _qualitieswin = _mm_insert_epi16(_qualitieswin, 64 << 2, 0);
        }
        
        // S even
        
        _m1 = _mm_or_si128(_initmask2, _mm_andnot_si128(_initmask, _m1));
        _m2 = _mm_or_si128(_initmask2, _mm_andnot_si128(_initmask, _m2));
        _m1 = _mm_min_epi16(_m1, _mm_min_epi16(_i1, _d1));
        
        if (s / 2 >= target_len) {
            minscore = std::min(static_cast<short>(extract_epi16(_m1, s / 2 - target_len)), minscore);
        }
        
        _m1 = _mm_add_epi16(_m1, _mm_min_epi16(_mm_andnot_si128(_mm_cmpeq_epi16(_targetwin, _truthwin),
                                                                _qualitieswin), _truthnqual));
        _d1 = _mm_min_epi16(_mm_add_epi16(_d2, _gap_extend),
                            _mm_add_epi16(_mm_min_epi16(_m2, _i2),
                                          _mm_srli_si128(_gap_open, 2))); // allow I->D
        _d1 = _mm_insert_epi16(_mm_slli_si128(_d1, 2), inf, 0);
        _i1 = _mm_add_epi16(_mm_min_epi16(_mm_add_epi16(_i2, _gap_extend),
                                          _mm_add_epi16(_m2, _gap_open)),
                            _nuc_prior);
        
        // S odd
        // truth needs updating; target is current
        const auto pos = bandSize + s / 2;
        const char base {(pos < truth_len) ? truth[pos] : 'N'};
        
        _truthwin   = _mm_insert_epi16(_mm_srli_si128(_truthwin, 2), base,
                                       bandSize - 1);
        _truthnqual = _mm_insert_epi16(_mm_srli_si128(_truthnqual, 2), base == 'N' ? nScore : inf,
                                       bandSize - 1);
        _gap_open   = _mm_insert_epi16(_mm_srli_si128(_gap_open, 2),
                                       gap_open[pos < truth_len ? pos : truth_len - 1] << 2,
                                       bandSize - 1);
        
        _initmask  = _mm_slli_si128(_initmask, 2);
        _initmask2 = _mm_slli_si128(_initmask2, 2);
        
        _m2 = _mm_min_epi16(_m2, _mm_min_epi16(_i2, _d2));
        
        if (s / 2 >= target_len) {
            minscore = std::min(static_cast<short>(extract_epi16(_m2, s / 2 - target_len)), minscore);
        }
        
        _m2 = _mm_add_epi16(_m2, _mm_min_epi16(_mm_andnot_si128(_mm_cmpeq_epi16(_targetwin, _truthwin),
                                                                 _qualitieswin), _truthnqual));
        _d2 = _mm_min_epi16(_mm_add_epi16(_d1, _gap_extend),
                            _mm_add_epi16(_mm_min_epi16(_m1, _i1), _gap_open)); // allow I->D
        _i2 = _mm_insert_epi16(_mm_add_epi16(_mm_min_epi16(_mm_add_epi16(_mm_srli_si128(_i1, 2),
                                                                         _gap_extend),
                                                           _mm_add_epi16(_mm_srli_si128(_m1, 2),
                                                                         _gap_open)),
                                             _nuc_prior), inf, bandSize - 1);
        
    }
    
    return (minscore + 0x8000) >> 2;
}

int align(const char* truth, const char* target, const std::int8_t* qualities,
          const int truth_len, const int target_len,
          const char* snv_mask, const std::int8_t* snv_prior,
          const std::int8_t* gap_open, short gap_extend, short nuc_prior) noexcept
{
    assert(truth_len > bandSize && (truth_len == target_len + 2 * bandSize - 1));
    
    gap_extend <<= 2;
    nuc_prior <<= 2;
    
    using SimdInt = __m128i;
    
    SimdInt _m1 {_mm_set1_epi16(inf)};
    auto _i1 = _m1;
    auto _d1 = _m1;
    auto _m2 = _m1;
    auto _i2 = _m1;
    auto _d2 = _m1;
    
    SimdInt _gap_extend {_mm_set1_epi16(gap_extend)};
    SimdInt _nuc_prior  {_mm_set1_epi16(nuc_prior)};
    SimdInt _initmask   {_mm_set_epi16(0,0,0,0,0,0,0,-1)};
    SimdInt _initmask2  {_mm_set_epi16(0,0,0,0,0,0,0,-0x8000)};
    
    SimdInt _truthwin  {_mm_set_epi16(truth[7], truth[6], truth[5], truth[4],
                                      truth[3], truth[2], truth[1], truth[0])};
    SimdInt _targetwin  {_m1};
    SimdInt _qualitieswin {_mm_set1_epi16(64 << 2)};
    
    SimdInt _snvmaskwin  {_mm_set_epi16(snv_mask[7], snv_mask[6], snv_mask[5], snv_mask[4],
                                        snv_mask[3], snv_mask[2], snv_mask[1], snv_mask[0])};
    SimdInt _snv_priorwin {_mm_set_epi16(snv_prior[7] << 2, snv_prior[6] << 2, snv_prior[5] << 2, snv_prior[4] << 2,
                                         snv_prior[3] << 2, snv_prior[2] << 2, snv_prior[1] << 2, snv_prior[0] << 2)};
    
    SimdInt _snvmask;
    
    // if N, make nScore; if != N, make inf
    SimdInt _truthnqual {_mm_add_epi16(_mm_and_si128(_mm_cmpeq_epi16(_truthwin, _mm_set1_epi16('N')),
                                                     _mm_set1_epi16(nScore - inf)),
                                       _mm_set1_epi16(inf))};
    
    SimdInt _gap_open {_mm_set_epi16(gap_open[7] << 2,gap_open[6] << 2,gap_open[5] << 2,gap_open[4] << 2,
                                     gap_open[3] << 2,gap_open[2] << 2,gap_open[1] << 2,gap_open[0] << 2)};
    
    short minscore {inf};
    
    for (int s {0}; s <= 2 * (target_len + bandSize); s += 2) {
        // truth is current; target needs updating
        _targetwin    = _mm_slli_si128(_targetwin, 2);
        _qualitieswin = _mm_slli_si128(_qualitieswin, 2);
        
        if (s / 2 < target_len) {
            _targetwin    = _mm_insert_epi16(_targetwin, target[s / 2], 0);
            _qualitieswin = _mm_insert_epi16(_qualitieswin, qualities[s / 2] << 2, 0);
        } else {
            _targetwin    = _mm_insert_epi16(_targetwin, '0', 0);
            _qualitieswin = _mm_insert_epi16(_qualitieswin, 64 << 2, 0);
        }
        
        // S even
        
        _m1 = _mm_or_si128(_initmask2, _mm_andnot_si128(_initmask, _m1));
        _m2 = _mm_or_si128(_initmask2, _mm_andnot_si128(_initmask, _m2));
        _m1 = _mm_min_epi16(_m1, _mm_min_epi16(_i1, _d1));
        
        if (s / 2 >= target_len) {
            minscore = std::min(static_cast<short>(extract_epi16(_m1, s / 2 - target_len)), minscore);
        }
        
        _snvmask = _mm_cmpeq_epi16(_targetwin, _snvmaskwin);
        
        _m1 = _mm_add_epi16(_m1, _mm_min_epi16(_mm_andnot_si128(_mm_cmpeq_epi16(_targetwin, _truthwin),
                                                                 _mm_min_epi16(_qualitieswin,
                                                                               _mm_or_si128(_mm_and_si128(_snvmask, _snv_priorwin),
                                                                                            _mm_andnot_si128(_snvmask, _qualitieswin)))),
                                               _truthnqual));
        _d1 = _mm_min_epi16(_mm_add_epi16(_d2, _gap_extend),
                            _mm_add_epi16(_mm_min_epi16(_m2, _i2),
                                          _mm_srli_si128(_gap_open, 2))); // allow I->D
        _d1 = _mm_insert_epi16(_mm_slli_si128(_d1, 2), inf, 0);
        _i1 = _mm_add_epi16(_mm_min_epi16(_mm_add_epi16(_i2, _gap_extend),
                                          _mm_add_epi16(_m2, _gap_open)),
                            _nuc_prior);
        
        // S odd
        // truth needs updating; target is current
        const auto pos = bandSize + s / 2;
        const char base {pos < truth_len ? truth[pos] : 'N'};
        
        _truthwin     = _mm_insert_epi16(_mm_srli_si128(_truthwin, 2),
                                         base,
                                         bandSize - 1);
        _truthnqual   = _mm_insert_epi16(_mm_srli_si128(_truthnqual, 2),
                                         base == 'N' ? nScore : inf,
                                         bandSize - 1);
        _snvmaskwin   = _mm_insert_epi16(_mm_srli_si128(_snvmaskwin, 2),
                                         pos < truth_len ? snv_mask[pos] : 'N',
                                         bandSize - 1);
        _snv_priorwin = _mm_insert_epi16(_mm_srli_si128(_snv_priorwin, 2),
                                         (pos < truth_len ? snv_prior[pos] : inf) << 2,
                                         bandSize - 1);
        _gap_open     = _mm_insert_epi16(_mm_srli_si128(_gap_open, 2),
                                         gap_open[pos < truth_len ? pos : truth_len - 1] << 2,
                                         bandSize - 1);
        
        _initmask  = _mm_slli_si128(_initmask, 2);
        _initmask2 = _mm_slli_si128(_initmask2, 2);
        
        _m2 = _mm_min_epi16(_m2, _mm_min_epi16(_i2, _d2));
        
        if (s / 2 >= target_len) {
            minscore = std::min(static_cast<short>(extract_epi16(_m2, s / 2 - target_len)), minscore);
        }
        
        _snvmask = _mm_cmpeq_epi16(_targetwin, _snvmaskwin);
        
        _m2 = _mm_add_epi16(_m2, _mm_min_epi16(_mm_andnot_si128(_mm_cmpeq_epi16(_targetwin, _truthwin),
                                                                _mm_min_epi16(_qualitieswin,
                                                                              _mm_or_si128(_mm_and_si128(_snvmask, _snv_priorwin),
                                                                                           _mm_andnot_si128(_snvmask, _qualitieswin)))),
                                               _truthnqual));
        _d2 = _mm_min_epi16(_mm_add_epi16(_d1, _gap_extend),
                            _mm_add_epi16(_mm_min_epi16(_m1, _i1), _gap_open)); // allow I->D
        _i2 = _mm_insert_epi16(_mm_add_epi16(_mm_min_epi16(_mm_add_epi16(_mm_srli_si128(_i1, 2),
                                                                         _gap_extend),
                                                           _mm_add_epi16(_mm_srli_si128(_m1, 2),
                                                                         _gap_open)),
                                             _nuc_prior), inf, bandSize - 1);
    }
    
    return (minscore + 0x8000) >> 2;
}

int align(const char* truth, const char* target, const std::int8_t* qualities,
          const int truth_len, const int target_len,
          const std::int8_t* gap_open, short gap_extend, short nuc_prior,
          int& first_pos, char* aln1, char* aln2) noexcept
{
    // target is the read; the shorter of the sequences
    // no checks for overflow are done
    
    // the bottom-left and top-right corners of the DP table are just
    // included at the extreme ends of the diagonal, which measures
    // n=8 entries diagonally across.  This fixes the length of the
    // longer (horizontal) sequence to 15 (2*8-1) more than the shorter
    
    assert(truth_len > bandSize && (truth_len == target_len + 2 * bandSize - 1));
    assert(aln1 != nullptr && aln2 != nullptr);
    
    gap_extend <<= 2;
    nuc_prior <<= 2;
    
    constexpr int matchLabel  {0};
    constexpr int insertLabel {1};
    constexpr int deleteLabel {3};
    
    using SimdInt = __m128i;
    
    SimdInt _m1 {_mm_set1_epi16(inf)};
    auto _i1 = _m1;
    auto _d1 = _m1;
    auto _m2 = _m1;
    auto _i2 = _m1;
    auto _d2 = _m1;
    
    SimdInt _gap_extend {_mm_set1_epi16(gap_extend)};
    SimdInt _nuc_prior  {_mm_set1_epi16(nuc_prior)};
    SimdInt _initmask   {_mm_set_epi16(0,0,0,0,0,0,0,-1)};
    SimdInt _initmask2  {_mm_set_epi16(0,0,0,0,0,0,0,-0x8000)};
    
    static const SimdInt _three {_mm_set1_epi16(3)};
    SmallVector<SimdInt> _backpointers(2 * (truth_len + bandSize));
    
    // sequence 1 is initialized with the n-long prefix, in forward direction
    // sequence 2 is initialized as empty; reverse direction
    SimdInt _truthwin  {_mm_set_epi16(truth[7], truth[6], truth[5], truth[4],
                                      truth[3], truth[2], truth[1], truth[0])};
    SimdInt _targetwin  {_m1};
    SimdInt _qualitieswin {_mm_set1_epi16(64 << 2)};
    
    // if N, make nScore; if != N, make inf
    SimdInt _truthnqual {_mm_add_epi16(_mm_and_si128(_mm_cmpeq_epi16(_truthwin, _mm_set1_epi16('N')),
                                                     _mm_set1_epi16(nScore - inf)),
                                       _mm_set1_epi16(inf))};
    
    SimdInt _gap_open {_mm_set_epi16(gap_open[7] << 2,gap_open[6] << 2,gap_open[5] << 2,gap_open[4] << 2,
                                     gap_open[3] << 2,gap_open[2] << 2,gap_open[1] << 2,gap_open[0] << 2)};
    
    short cur_score {0}, minscore {inf}, minscoreidx {-1};
    
    // main loop.  Do one extra iteration, with nucs from sequence 2 just moved out
    // of the targetwin/qual arrays, to simplify getting back pointers
    int s;
    for (s = 0; s <= 2 * (target_len + bandSize); s += 2) {
        // truth is current; target needs updating
        _targetwin    = _mm_slli_si128(_targetwin, 2);
        _qualitieswin = _mm_slli_si128(_qualitieswin, 2);
        
        if (s / 2 < target_len) {
            _targetwin    = _mm_insert_epi16(_targetwin, target[s / 2], 0);
            _qualitieswin = _mm_insert_epi16(_qualitieswin, qualities[s / 2] << 2, 0);
        } else {
            _targetwin    = _mm_insert_epi16(_targetwin, '0', 0);
            _qualitieswin = _mm_insert_epi16(_qualitieswin, 64 << 2, 0);
        }
        
        // S even
        _m1 = _mm_or_si128(_initmask2, _mm_andnot_si128(_initmask, _m1));
        _m2 = _mm_or_si128(_initmask2, _mm_andnot_si128(_initmask, _m2));
        _m1 = _mm_min_epi16(_m1, _mm_min_epi16(_i1, _d1));
        
        if (s / 2 >= target_len) {
            cur_score = extract_epi16(_m1, s / 2 - target_len);
            if (cur_score < minscore) {
                minscore = cur_score;
                minscoreidx = s;     // point back to the match state at this entry, so as not to
            }                        // have to store the state at s-2
        }
        
        _m1 = _mm_add_epi16(_m1, _mm_min_epi16(_mm_andnot_si128(_mm_cmpeq_epi16(_targetwin, _truthwin),
                                                                _qualitieswin), _truthnqual));
        _d1 = _mm_min_epi16(_mm_add_epi16(_d2, _gap_extend),
                            _mm_add_epi16(_mm_min_epi16(_m2, _i2),
                                          _mm_srli_si128(_gap_open, 2))); // allow I->D
        _d1 = _mm_insert_epi16(_mm_slli_si128(_d1, 2), inf, 0);
        _i1 = _mm_add_epi16(_mm_min_epi16(_mm_add_epi16(_i2, _gap_extend),
                                          _mm_add_epi16(_m2, _gap_open)), _nuc_prior);
        
        _backpointers[s] = _mm_or_si128(_mm_or_si128(_mm_and_si128(_three, _m1),
                                                     _mm_slli_epi16(_mm_and_si128(_three, _i1), 2 * insertLabel)),
                                        _mm_slli_epi16(_mm_and_si128(_three, _d1), 2 * deleteLabel));
        
        // set state labels
        _m1 = _mm_andnot_si128(_three, _m1);
        _i1 = _mm_or_si128(_mm_andnot_si128(_three, _i1), _mm_srli_epi16(_three, 1));
        _d1 = _mm_or_si128(_mm_andnot_si128(_three, _d1), _three);
        
        // S odd
        
        // truth needs updating; target is current
        const char c {(bandSize + s / 2 < truth_len) ? truth[bandSize + (s / 2)] : 'N'};
        
        _truthwin   = _mm_insert_epi16(_mm_srli_si128(_truthwin,   2), c, bandSize - 1);
        _truthnqual = _mm_insert_epi16(_mm_srli_si128(_truthnqual, 2), (c == 'N') ? nScore : inf, bandSize - 1);
        _gap_open   = _mm_insert_epi16(_mm_srli_si128(_gap_open,  2),
                                       gap_open[bandSize + s / 2 < truth_len ? bandSize + s / 2 : truth_len - 1] << 2, bandSize - 1);
        _initmask  = _mm_slli_si128(_initmask, 2);
        _initmask2 = _mm_slli_si128(_initmask2, 2);
        _m2 = _mm_min_epi16(_m2, _mm_min_epi16(_i2, _d2));
        
        // at this point, extract minimum score.  Referred-to position must
        // be y==target_len-1, so that current position has y==target_len; i==0 so d=0 and y=s/2
        if (s / 2 >= target_len) {
            cur_score = extract_epi16(_m2, s / 2 - target_len);
            if (cur_score < minscore) {
                minscore = cur_score;
                minscoreidx = s + 1;
            }
        }
        
        _m2 = _mm_add_epi16(_m2, _mm_min_epi16(_mm_andnot_si128(_mm_cmpeq_epi16(_targetwin, _truthwin), _qualitieswin),
                                               _truthnqual));
        _d2 = _mm_min_epi16(_mm_add_epi16(_d1, _gap_extend),
                            _mm_add_epi16(_mm_min_epi16(_m1, _i1),  // allow I->D
                                          _gap_open));
        _i2 = _mm_insert_epi16(_mm_add_epi16(_mm_min_epi16(_mm_add_epi16(_mm_srli_si128(_i1, 2), _gap_extend),
                                                           _mm_add_epi16(_mm_srli_si128(_m1, 2), _gap_open)),
                                             _nuc_prior),
                               inf, bandSize - 1);
        _backpointers[s + 1] = _mm_or_si128(_mm_or_si128(_mm_and_si128(_three, _m2),
                                                         _mm_slli_epi16(_mm_and_si128(_three, _i2), 2 * insertLabel)),
                                            _mm_slli_epi16(_mm_and_si128(_three, _d2), 2 * deleteLabel));
        
        // set state labels
        _m2 = _mm_andnot_si128(_three, _m2);
        _i2 = _mm_or_si128(_mm_andnot_si128(_three, _i2), _mm_srli_epi16(_three, 1));
        _d2 = _mm_or_si128(_mm_andnot_si128(_three, _d2), _three);
    }
    
    s = minscoreidx;    // point to the dummy match transition
    
    auto i      = s / 2 - target_len;
    auto y      = target_len;
    auto x      = s - y;
    auto alnidx = 0;
    auto state = (reinterpret_cast<short*>(_backpointers.data() + s)[i] >> (2 * matchLabel)) & 3;
    
    s -= 2;
    
    // this is 2*y (s even) or 2*y+1 (s odd)
    while (y > 0) {
        const auto new_state = (reinterpret_cast<short*>(_backpointers.data() + s)[i] >> (2 * state)) & 3;
        
        if (state == matchLabel) {
            s -= 2;
            aln1[alnidx] = truth[--x];
            aln2[alnidx] = target[--y];
        } else if (state == insertLabel) {
            i += s & 1;
            s -= 1;
            aln1[alnidx] = gap;
            aln2[alnidx] = target[--y];
        } else {
            s -= 1;
            i -= s & 1;
            aln1[alnidx] = truth[--x];
            aln2[alnidx] = gap;
        }
        state = new_state;
        alnidx++;
    }
    
    aln1[alnidx] = 0;
    aln2[alnidx] = 0;
    
    first_pos = x;
    
    // reverse them
    for (int j {alnidx - 1}, i = 0; i < j; ++i, j--) {
        x = aln1[i];
        y = aln2[i];
        aln1[i] = aln1[j];
        aln2[i] = aln2[j];
        aln1[j] = x;
        aln2[j] = y;
    }
    
    return (minscore + 0x8000) >> 2;
}

int align(const char* truth, const char* target, const std::int8_t* qualities,
          int truth_len, int target_len,
          const char* snv_mask, const std::int8_t* snv_prior,
          const std::int8_t* gap_open, short gap_extend, short nuc_prior,
          char* aln1, char* aln2, int& first_pos) noexcept
{
    assert(truth_len > bandSize && (truth_len == target_len + 2 * bandSize - 1));
    assert(aln1 != nullptr && aln2 != nullptr);
    
    gap_extend <<= 2;
    nuc_prior <<= 2;
    
    constexpr int matchLabel  {0};
    constexpr int insertLabel {1};
    constexpr int deleteLabel {3};
    
    using SimdInt = __m128i;
    
    SimdInt _m1 {_mm_set1_epi16(inf)};
    auto _i1 = _m1;
    auto _d1 = _m1;
    auto _m2 = _m1;
    auto _i2 = _m1;
    auto _d2 = _m1;
    
    SimdInt _gap_extend {_mm_set1_epi16(gap_extend)};
    SimdInt _nuc_prior  {_mm_set1_epi16(nuc_prior)};
    SimdInt _initmask   {_mm_set_epi16(0,0,0,0,0,0,0,-1)};
    SimdInt _initmask2  {_mm_set_epi16(0,0,0,0,0,0,0,-0x8000)};
    
    static const SimdInt _three {_mm_set1_epi16(3)};
    SmallVector<SimdInt> _backpointers(2 * (truth_len + bandSize));
    
    SimdInt _truthwin  {_mm_set_epi16(truth[7], truth[6], truth[5], truth[4],
                                      truth[3], truth[2], truth[1], truth[0])};
    SimdInt _targetwin  {_m1};
    SimdInt _qualitieswin {_mm_set1_epi16(64 << 2)};
    
    SimdInt _snvmaskwin  {_mm_set_epi16(snv_mask[7], snv_mask[6], snv_mask[5], snv_mask[4],
                                        snv_mask[3], snv_mask[2], snv_mask[1], snv_mask[0])};
    SimdInt _snv_priorwin {_mm_set_epi16(snv_prior[7] << 2, snv_prior[6] << 2, snv_prior[5] << 2, snv_prior[4] << 2,
                                         snv_prior[3] << 2, snv_prior[2] << 2, snv_prior[1] << 2, snv_prior[0] << 2)};
    
    SimdInt _snvmask;
    
    // if N, make nScore; if != N, make inf
    SimdInt _truthnqual {_mm_add_epi16(_mm_and_si128(_mm_cmpeq_epi16(_truthwin, _mm_set1_epi16('N')),
                                                     _mm_set1_epi16(nScore - inf)),
                                       _mm_set1_epi16(inf))};
    SimdInt _gap_open {_mm_set_epi16(gap_open[7] << 2,gap_open[6] << 2,gap_open[5] << 2,gap_open[4] << 2,
                                     gap_open[3] << 2,gap_open[2] << 2,gap_open[1] << 2,gap_open[0] << 2)};
    
    short cur_score {0}, minscore {inf}, minscoreidx {-1};
    
    int s;
    for (s = 0; s <= 2 * (target_len + bandSize); s += 2) {
        // truth is current; target needs updating
        _targetwin    = _mm_slli_si128(_targetwin, 2);
        _qualitieswin = _mm_slli_si128(_qualitieswin, 2);
        
        if (s / 2 < target_len) {
            _targetwin    = _mm_insert_epi16(_targetwin, target[s / 2], 0);
            _qualitieswin = _mm_insert_epi16(_qualitieswin, qualities[s / 2] << 2, 0);
        } else {
            _targetwin    = _mm_insert_epi16(_targetwin, '0', 0);
            _qualitieswin = _mm_insert_epi16(_qualitieswin, 64 << 2, 0);
        }
        
        // S even
        
        // initialize to -0x8000
        _m1 = _mm_or_si128(_initmask2, _mm_andnot_si128(_initmask, _m1));
        _m2 = _mm_or_si128(_initmask2, _mm_andnot_si128(_initmask, _m2));
        _m1 = _mm_min_epi16(_m1, _mm_min_epi16(_i1, _d1));
        
        // at this point, extract minimum score.  Referred-to position must
        // be y==target_len-1, so that current position has y==target_len; i==0 so d=0 and y=s/2
        
        if (s / 2 >= target_len) {
            cur_score = extract_epi16(_m1, s / 2 - target_len);
            if (cur_score < minscore) {
                minscore = cur_score;
                minscoreidx = s;     // point back to the match state at this entry, so as not to
            }                        // have to store the state at s-2
        }
        
        _snvmask = _mm_cmpeq_epi16(_targetwin, _snvmaskwin);
        _m1 = _mm_add_epi16(_m1, _mm_min_epi16(_mm_andnot_si128(_mm_cmpeq_epi16(_targetwin, _truthwin),
                                                                _mm_min_epi16(_qualitieswin,
                                                                              _mm_or_si128(_mm_and_si128(_snvmask, _snv_priorwin),
                                                                                           _mm_andnot_si128(_snvmask, _qualitieswin)))),
                                               _truthnqual));
        _d1 = _mm_min_epi16(_mm_add_epi16(_d2, _gap_extend),
                            _mm_add_epi16(_mm_min_epi16(_m2, _i2),
                                          _mm_srli_si128(_gap_open, 2))); // allow I->D
        _d1 = _mm_insert_epi16(_mm_slli_si128(_d1, 2), inf, 0);
        _i1 = _mm_add_epi16(_mm_min_epi16(_mm_add_epi16(_i2, _gap_extend),
                                          _mm_add_epi16(_m2, _gap_open)), _nuc_prior);
        
        _backpointers[s] = _mm_or_si128(_mm_or_si128(_mm_and_si128(_three, _m1),
                                                     _mm_slli_epi16(_mm_and_si128(_three, _i1),
                                                                    2 * insertLabel)),
                                        _mm_slli_epi16(_mm_and_si128(_three, _d1), 2 * deleteLabel));
        
        // set state labels
        _m1 = _mm_andnot_si128(_three, _m1);
        _i1 = _mm_or_si128(_mm_andnot_si128(_three, _i1), _mm_srli_epi16(_three, 1));
        _d1 = _mm_or_si128(_mm_andnot_si128(_three, _d1), _three);
        
        // S odd
        
        // truth needs updating; target is current
        const auto pos = bandSize + s / 2;
        const char base {(pos < truth_len) ? truth[pos] : 'N'};
        
        _truthwin   = _mm_insert_epi16(_mm_srli_si128(_truthwin, 2), base, bandSize - 1);
        _truthnqual = _mm_insert_epi16(_mm_srli_si128(_truthnqual, 2),
                                       (base == 'N') ? nScore : inf, bandSize - 1);
        _snvmaskwin   = _mm_insert_epi16(_mm_srli_si128(_snvmaskwin, 2),
                                         pos < truth_len ? snv_mask[pos] : 'N', bandSize - 1);
        _snv_priorwin = _mm_insert_epi16(_mm_srli_si128(_snv_priorwin, 2),
                                         (pos < truth_len) ? snv_prior[pos] << 2 : inf << 2,
                                         bandSize - 1);
        _gap_open  = _mm_insert_epi16(_mm_srli_si128(_gap_open,  2),
                                      gap_open[pos < truth_len ? pos : truth_len - 1] << 2,
                                      bandSize - 1);
        
        _initmask  = _mm_slli_si128(_initmask, 2);
        _initmask2 = _mm_slli_si128(_initmask2, 2);
        
        _m2 = _mm_min_epi16(_m2, _mm_min_epi16(_i2, _d2));
        
        // at this point, extract minimum score.  Referred-to position must
        // be y==target_len-1, so that current position has y==target_len; i==0 so d=0 and y=s/2
        if (s / 2 >= target_len) {
            cur_score = extract_epi16(_m2, s / 2 - target_len);
            if (cur_score < minscore) {
                minscore = cur_score;
                minscoreidx = s + 1;
            }
        }
        
        _snvmask = _mm_cmpeq_epi16(_targetwin, _snvmaskwin);
        
        _m2 = _mm_add_epi16(_m2, _mm_min_epi16(_mm_andnot_si128(_mm_cmpeq_epi16(_targetwin, _truthwin),
                                                                _mm_min_epi16(_qualitieswin,
                                                                              _mm_or_si128(_mm_and_si128(_snvmask, _snv_priorwin),
                                                                                           _mm_andnot_si128(_snvmask, _qualitieswin)))),
                                               _truthnqual));
        _d2 = _mm_min_epi16(_mm_add_epi16(_d1, _gap_extend),
                            _mm_add_epi16(_mm_min_epi16(_m1, _i1),  // allow I->D
                                          _gap_open));
        _i2 = _mm_insert_epi16(_mm_add_epi16(_mm_min_epi16(_mm_add_epi16(_mm_srli_si128(_i1, 2), _gap_extend),
                                                           _mm_add_epi16(_mm_srli_si128(_m1, 2), _gap_open)),
                                             _nuc_prior),
                               inf, bandSize - 1);
        _backpointers[s + 1] = _mm_or_si128(_mm_or_si128(_mm_and_si128(_three, _m2),
                                                         _mm_slli_epi16(_mm_and_si128(_three, _i2), 2 * insertLabel)),
                                            _mm_slli_epi16(_mm_and_si128(_three, _d2), 2 * deleteLabel));
        
        // set state labels
        _m2 = _mm_andnot_si128(_three, _m2);
        _i2 = _mm_or_si128(_mm_andnot_si128(_three, _i2), _mm_srli_epi16(_three, 1));
        _d2 = _mm_or_si128(_mm_andnot_si128(_three, _d2), _three);
    }
    
    s = minscoreidx;    // point to the dummy match transition
    
    auto i      = s / 2 - target_len;
    auto y      = target_len;
    auto x      = s - y;
    auto alnidx = 0;
    auto state  = (reinterpret_cast<short*>(_backpointers.data() + s)[i] >> (2 * matchLabel)) & 3;
    
    s -= 2;
    
    // this is 2*y (s even) or 2*y+1 (s odd)
    while (y > 0) {
        const auto new_state = (reinterpret_cast<short*>(_backpointers.data() + s)[i] >> (2 * state)) & 3;
        
        if (state == matchLabel) {
            s -= 2;
            aln1[alnidx] = truth[--x];
            aln2[alnidx] = target[--y];
        } else if (state == insertLabel) {
            i += s & 1;
            s -= 1;
            aln1[alnidx] = gap;
            aln2[alnidx] = target[--y];
        } else {
            s -= 1;
            i -= s & 1;
            aln1[alnidx] = truth[--x];
            aln2[alnidx] = gap;
        }
        state = new_state;
        alnidx++;
    }
    
    aln1[alnidx] = 0;
    aln2[alnidx] = 0;
    first_pos = x;
    
    // reverse them
    for (int j {alnidx - 1}, i = 0; i < j; ++i, j--) {
        x = aln1[i];
        y = aln2[i];
        aln1[i] = aln1[j];
        aln2[i] = aln2[j];
        aln1[j] = x;
        aln2[j] = y;
    }
    
    return (minscore + 0x8000) >> 2;
}

int calculate_flank_score(const int truth_len, const int lhs_flank_len, const int rhs_flank_len,
                          const std::int8_t* quals, const std::int8_t* gap_open,
                          const short gap_extend, const short nuc_prior,
                          const int first_pos, const char* aln1, const char* aln2,
                          int& target_mask_size) noexcept
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
        if (aln1[i] == gap) {
            new_state = insertion;
        } else if (aln2[i] == gap) { // can't be both '-'
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
                            result += nScore >> 2;
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

int calculate_flank_score(const int truth_len, const int lhs_flank_len, const int rhs_flank_len,
                          const char* target, const std::int8_t* quals,
                          const char* snv_mask, const std::int8_t* snv_prior,
                          const std::int8_t* gap_open, const short gap_extend, const short nuc_prior,
                          const int first_pos, const char* aln1, const char* aln2,
                          int& target_mask_size) noexcept
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
        if (aln1[i] == gap) {
            new_state = insertion;
        } else if (aln2[i] == gap) { // can't be both '-'
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
                            result += nScore >> 2;
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

} // namespace simd
} // namespace hmm
} // namespace octopus
