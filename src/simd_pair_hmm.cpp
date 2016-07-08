/*****************************************************************************************************************
 This code is copyright (c) Gerton Lunter, Dec 2009
 It may not be distributed, made public, or used in other software without the permission of the copyright holder
 ******************************************************************************************************************/

#include "simd_pair_hmm.hpp"

#include <emmintrin.h>
#include <algorithm>
#include <cassert>

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

namespace SimdPairHmm
{
constexpr short N_SCORE {0 << 2};
constexpr int BAND_SIZE {8};
constexpr short INF {0x7800};

constexpr char Gap {'-'};

int align(const char* truth, const char* target, const std::int8_t* qualities,
          const int truth_len, const int target_len,
          const std::int8_t* gap_open, short gap_extend, short nuc_prior)
{
    // target is the read; the shorter of the sequences
    // no checks for overflow are done
    
    // the bottom-left and top-right corners of the DP table are just
    // included at the extreme ends of the diagonal, which measures
    // n=8 entries diagonally across.  This fixes the length of the
    // longer (horizontal) sequence to 15 (2*8-1) more than the shorter
    
    // the << 2's are because the lower two bits are reserved for back tracing
    
    assert(truth_len > BAND_SIZE && (truth_len == target_len + 2 * BAND_SIZE - 1));
    
    gap_extend <<= 2;
    nuc_prior <<= 2;
    
    using SimdInt = __m128i;
    
    SimdInt _m1 {_mm_set1_epi16(INF)};
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
    
    // if N, make N_SCORE; if != N, make INF
    SimdInt _truthnqual {_mm_add_epi16(_mm_and_si128(_mm_cmpeq_epi16(_truthwin, _mm_set1_epi16('N')),
                                                     _mm_set1_epi16(N_SCORE - INF)),
                                       _mm_set1_epi16(INF))};
    
    SimdInt _gap_open {_mm_set_epi16(gap_open[7] << 2,gap_open[6] << 2,gap_open[5] << 2,gap_open[4] << 2,
                                     gap_open[3] << 2,gap_open[2] << 2,gap_open[1] << 2,gap_open[0] << 2)};
    
    short cur_score {0}, minscore {INF};
    
    for (int s {0}; s <= 2 * (target_len + BAND_SIZE); s += 2) {
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
        
        // at this point, extract minimum score.  Referred-to position must
        // be y==target_len-1, so that current position has y==target_len; i==0 so d=0 and y=s/2
        
        if (s / 2 >= target_len) {
            cur_score = _mm_extract_epi16(_m1, s / 2 - target_len);
            
            if (cur_score < minscore) {
                minscore = cur_score;
            }
        }
        
        _m1 = _mm_add_epi16(_m1, _mm_min_epi16(_mm_andnot_si128(_mm_cmpeq_epi16(_targetwin, _truthwin),
                                                                _qualitieswin), _truthnqual));
        
        _d1 = _mm_min_epi16(_mm_add_epi16(_d2, _gap_extend),
                            _mm_add_epi16(_mm_min_epi16(_m2, _i2),
                                          _mm_srli_si128(_gap_open, 2))); // allow I->D
        
        _d1 = _mm_insert_epi16(_mm_slli_si128(_d1, 2), INF, 0);
        
        _i1 = _mm_add_epi16(_mm_min_epi16(_mm_add_epi16(_i2, _gap_extend),
                                          _mm_add_epi16(_m2, _gap_open)),
                            _nuc_prior);
        
        // S odd
        
        // truth needs updating; target is current
        const auto pos = BAND_SIZE + s / 2;
        
        const char base {(pos < truth_len) ? truth[pos] : 'N'};
        
        _truthwin   = _mm_insert_epi16(_mm_srli_si128(_truthwin, 2), base,
                                       BAND_SIZE - 1);
        _truthnqual = _mm_insert_epi16(_mm_srli_si128(_truthnqual, 2), base == 'N' ? N_SCORE : INF,
                                       BAND_SIZE - 1);
        _gap_open   = _mm_insert_epi16(_mm_srli_si128(_gap_open, 2),
                                       gap_open[pos < truth_len ? pos : truth_len - 1] << 2,
                                       BAND_SIZE - 1);
        
        _initmask  = _mm_slli_si128(_initmask, 2);
        _initmask2 = _mm_slli_si128(_initmask2, 2);
        
        _m2 = _mm_min_epi16(_m2, _mm_min_epi16(_i2, _d2));
        
        // at this point, extract minimum score.  Referred-to position must
        // be y==target_len-1, so that current position has y==target_len; i==0 so d=0 and y=s/2
        if (s / 2 >= target_len) {
            cur_score = _mm_extract_epi16(_m2, s / 2 - target_len);
            
            if (cur_score < minscore) {
                minscore = cur_score;
            }
        }
        
        _m2 = _mm_add_epi16(_m2, _mm_min_epi16(_mm_andnot_si128(_mm_cmpeq_epi16(_targetwin, _truthwin),
                                                                 _qualitieswin), _truthnqual));
        
        _d2 = _mm_min_epi16(_mm_add_epi16(_d1, _gap_extend),
                            _mm_add_epi16(_mm_min_epi16(_m1, _i1), _gap_open)); // allow I->D
        
        _i2 = _mm_insert_epi16(_mm_add_epi16(_mm_min_epi16(_mm_add_epi16(_mm_srli_si128(_i1, 2),
                                                                         _gap_extend),
                                                           _mm_add_epi16(_mm_srli_si128(_m1, 2),
                                                                         _gap_open)),
                                             _nuc_prior), INF, BAND_SIZE - 1);
        
    }
    
    return (minscore + 0x8000) >> 2;
}

int align(const char* truth, const char* target, const std::int8_t* qualities,
          const int truth_len, const int target_len,
          const char* snv_mask, const std::int8_t* snv_prior,
          const std::int8_t* gap_open, short gap_extend, short nuc_prior)
{
    assert(truth_len > BAND_SIZE && (truth_len == target_len + 2 * BAND_SIZE - 1));
    
    gap_extend <<= 2;
    nuc_prior <<= 2;
    
    using SimdInt = __m128i;
    
    SimdInt _m1 {_mm_set1_epi16(INF)};
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
    
    // if N, make N_SCORE; if != N, make INF
    SimdInt _truthnqual {_mm_add_epi16(_mm_and_si128(_mm_cmpeq_epi16(_truthwin, _mm_set1_epi16('N')),
                                                     _mm_set1_epi16(N_SCORE - INF)),
                                       _mm_set1_epi16(INF))};
    
    SimdInt _gap_open {_mm_set_epi16(gap_open[7] << 2,gap_open[6] << 2,gap_open[5] << 2,gap_open[4] << 2,
                                     gap_open[3] << 2,gap_open[2] << 2,gap_open[1] << 2,gap_open[0] << 2)};
    
    short cur_score {0}, minscore {INF};
    
    for (int s {0}; s <= 2 * (target_len + BAND_SIZE); s += 2) {
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
            cur_score = _mm_extract_epi16(_m1, s / 2 - target_len);
            
            if (cur_score < minscore) {
                minscore = cur_score;
            }
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
        
        _d1 = _mm_insert_epi16(_mm_slli_si128(_d1, 2), INF, 0);
        
        _i1 = _mm_add_epi16(_mm_min_epi16(_mm_add_epi16(_i2, _gap_extend),
                                          _mm_add_epi16(_m2, _gap_open)),
                            _nuc_prior);
        
        // S odd
        
        // truth needs updating; target is current
        const auto pos = BAND_SIZE + s / 2;
        
        const char base {pos < truth_len ? truth[pos] : 'N'};
        
        _truthwin     = _mm_insert_epi16(_mm_srli_si128(_truthwin, 2),
                                         base,
                                         BAND_SIZE - 1);
        
        _truthnqual   = _mm_insert_epi16(_mm_srli_si128(_truthnqual, 2),
                                         base == 'N' ? N_SCORE : INF,
                                         BAND_SIZE - 1);
        
        _snvmaskwin   = _mm_insert_epi16(_mm_srli_si128(_snvmaskwin, 2),
                                         pos < truth_len ? snv_mask[pos] : 'N',
                                         BAND_SIZE - 1);
        
        _snv_priorwin = _mm_insert_epi16(_mm_srli_si128(_snv_priorwin, 2),
                                         (pos < truth_len ? snv_prior[pos] : INF) << 2,
                                         BAND_SIZE - 1);
        
        _gap_open     = _mm_insert_epi16(_mm_srli_si128(_gap_open, 2),
                                         gap_open[pos < truth_len ? pos : truth_len - 1] << 2,
                                         BAND_SIZE - 1);
        
        _initmask  = _mm_slli_si128(_initmask, 2);
        _initmask2 = _mm_slli_si128(_initmask2, 2);
        
        _m2 = _mm_min_epi16(_m2, _mm_min_epi16(_i2, _d2));
        
        if (s / 2 >= target_len) {
            cur_score = _mm_extract_epi16(_m2, s / 2 - target_len);
            
            if (cur_score < minscore) {
                minscore = cur_score;
            }
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
                                             _nuc_prior), INF, BAND_SIZE - 1);
        
    }
    
    return (minscore + 0x8000) >> 2;
}

int align(const char* truth, const char* target, const std::int8_t* qualities,
          const int truth_len, const int target_len,
          const std::int8_t* gap_open, short gap_extend, short nuc_prior,
          char* aln1, char* aln2, int& first_pos)
{
    // target is the read; the shorter of the sequences
    // no checks for overflow are done
    
    // the bottom-left and top-right corners of the DP table are just
    // included at the extreme ends of the diagonal, which measures
    // n=8 entries diagonally across.  This fixes the length of the
    // longer (horizontal) sequence to 15 (2*8-1) more than the shorter
    
    assert(truth_len > BAND_SIZE && (truth_len == target_len + 2 * BAND_SIZE - 1));
    
    assert(aln1 != nullptr && aln2 != nullptr);
    
    gap_extend <<= 2;
    nuc_prior <<= 2;
    
    constexpr int match_label  {0};
    constexpr int insert_label {1};
    constexpr int delete_label {3};
    
    using SimdInt = __m128i;
    
    SimdInt _m1 {_mm_set1_epi16(INF)};
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
    SimdInt _backpointers[2 * (truth_len + BAND_SIZE)];
    
    // sequence 1 is initialized with the n-long prefix, in forward direction
    // sequence 2 is initialized as empty; reverse direction
    SimdInt _truthwin  {_mm_set_epi16(truth[7], truth[6], truth[5], truth[4],
                                      truth[3], truth[2], truth[1], truth[0])};
    SimdInt _targetwin  {_m1};
    SimdInt _qualitieswin {_mm_set1_epi16(64 << 2)};
    
    // if N, make N_SCORE; if != N, make INF
    SimdInt _truthnqual {_mm_add_epi16(_mm_and_si128(_mm_cmpeq_epi16(_truthwin, _mm_set1_epi16('N')),
                                                     _mm_set1_epi16(N_SCORE - INF)),
                                       _mm_set1_epi16(INF))};
    
    SimdInt _gap_open {_mm_set_epi16(gap_open[7] << 2,gap_open[6] << 2,gap_open[5] << 2,gap_open[4] << 2,
                                     gap_open[3] << 2,gap_open[2] << 2,gap_open[1] << 2,gap_open[0] << 2)};
    
    short cur_score {0}, minscore {INF}, minscoreidx {-1};
    
    // main loop.  Do one extra iteration, with nucs from sequence 2 just moved out
    // of the targetwin/qual arrays, to simplify getting back pointers
    int s;
    
    for (s = 0; s <= 2 * (target_len + BAND_SIZE); s += 2) {
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
            cur_score = _mm_extract_epi16(_m1, s / 2 - target_len);
            
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
        
        _d1 = _mm_insert_epi16(_mm_slli_si128(_d1, 2), INF, 0);
        
        _i1 = _mm_add_epi16(_mm_min_epi16(_mm_add_epi16(_i2, _gap_extend),
                                          _mm_add_epi16(_m2, _gap_open)), _nuc_prior);
        
        _backpointers[s] = _mm_or_si128(_mm_or_si128(_mm_and_si128(_three, _m1),
                                                     _mm_slli_epi16(_mm_and_si128(_three, _i1), 2 * insert_label)),
                                        _mm_slli_epi16(_mm_and_si128(_three, _d1), 2 * delete_label));
        
        // set state labels
        _m1 = _mm_andnot_si128(_three, _m1);
        _i1 = _mm_or_si128(_mm_andnot_si128(_three, _i1), _mm_srli_epi16(_three, 1));
        _d1 = _mm_or_si128(_mm_andnot_si128(_three, _d1), _three);
        
        // S odd
        
        // truth needs updating; target is current
        const char c {(BAND_SIZE + s / 2 < truth_len) ? truth[BAND_SIZE + (s / 2)] : 'N'};
        
        _truthwin   = _mm_insert_epi16(_mm_srli_si128(_truthwin,   2), c, BAND_SIZE - 1);
        _truthnqual = _mm_insert_epi16(_mm_srli_si128(_truthnqual, 2), (c == 'N') ? N_SCORE : INF, BAND_SIZE - 1);
        _gap_open   = _mm_insert_epi16(_mm_srli_si128(_gap_open,  2),
                                       gap_open[BAND_SIZE + s / 2 < truth_len ? BAND_SIZE + s / 2 : truth_len - 1] << 2, BAND_SIZE - 1);
        
        _initmask  = _mm_slli_si128(_initmask, 2);
        _initmask2 = _mm_slli_si128(_initmask2, 2);
        
        _m2 = _mm_min_epi16(_m2, _mm_min_epi16(_i2, _d2));
        
        // at this point, extract minimum score.  Referred-to position must
        // be y==target_len-1, so that current position has y==target_len; i==0 so d=0 and y=s/2
        if (s / 2 >= target_len) {
            cur_score = _mm_extract_epi16(_m2, s / 2 - target_len);
            
            if (cur_score < minscore) {
                minscore = cur_score;
                minscoreidx = s+1;
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
                               INF, BAND_SIZE - 1);
        
        _backpointers[s + 1] = _mm_or_si128(_mm_or_si128(_mm_and_si128(_three, _m2),
                                                         _mm_slli_epi16(_mm_and_si128(_three, _i2), 2 * insert_label)),
                                            _mm_slli_epi16(_mm_and_si128(_three, _d2), 2 * delete_label));
        
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
    auto state  = ((((short*)(_backpointers + s))[i]) >> (2 * match_label)) & 3;
    
    s -= 2;
    
    // this is 2*y (s even) or 2*y+1 (s odd)
    while (y > 0) {
        const auto new_state = ((((short*)(_backpointers + s))[i]) >> (2 * state)) & 3;
        
        if (state == match_label) {
            s -= 2;
            aln1[alnidx] = truth[--x];
            aln2[alnidx] = target[--y];
        } else if (state == insert_label) {
            i += s & 1;
            s -= 1;
            aln1[alnidx] = Gap;
            aln2[alnidx] = target[--y];
        } else {
            s -= 1;
            i -= s & 1;
            aln1[alnidx] = truth[--x];
            aln2[alnidx] = Gap;
        }
        
        state = new_state;
        alnidx++;
    }
    
    aln1[alnidx] = 0;
    aln2[alnidx] = 0;
    
    first_pos = x;
    
    // reverse them
    int j;
    for (i = 0, j = alnidx - 1; i < j; ++i, j--) {
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
          char* aln1, char* aln2, int& first_pos)
{
    assert(truth_len > BAND_SIZE && (truth_len == target_len + 2 * BAND_SIZE - 1));
    assert(aln1 != nullptr && aln2 != nullptr);
    
    gap_extend <<= 2;
    nuc_prior <<= 2;
    
    constexpr int match_label  {0};
    constexpr int insert_label {1};
    constexpr int delete_label {3};
    
    using SimdInt = __m128i;
    
    SimdInt _m1 {_mm_set1_epi16(INF)};
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
    SimdInt _backpointers[2 * (truth_len + BAND_SIZE)];
    
    SimdInt _truthwin  {_mm_set_epi16(truth[7], truth[6], truth[5], truth[4],
                                      truth[3], truth[2], truth[1], truth[0])};
    SimdInt _targetwin  {_m1};
    SimdInt _qualitieswin {_mm_set1_epi16(64 << 2)};
    
    SimdInt _snvmaskwin  {_mm_set_epi16(snv_mask[7], snv_mask[6], snv_mask[5], snv_mask[4],
                                        snv_mask[3], snv_mask[2], snv_mask[1], snv_mask[0])};
    SimdInt _snv_priorwin {_mm_set_epi16(snv_prior[7] << 2, snv_prior[6] << 2, snv_prior[5] << 2, snv_prior[4] << 2,
                                         snv_prior[3] << 2, snv_prior[2] << 2, snv_prior[1] << 2, snv_prior[0] << 2)};
    
    SimdInt _snvmask;
    
    // if N, make N_SCORE; if != N, make INF
    SimdInt _truthnqual {_mm_add_epi16(_mm_and_si128(_mm_cmpeq_epi16(_truthwin, _mm_set1_epi16('N')),
                                                     _mm_set1_epi16(N_SCORE - INF)),
                                       _mm_set1_epi16(INF))};
    
    SimdInt _gap_open {_mm_set_epi16(gap_open[7] << 2,gap_open[6] << 2,gap_open[5] << 2,gap_open[4] << 2,
                                     gap_open[3] << 2,gap_open[2] << 2,gap_open[1] << 2,gap_open[0] << 2)};
    
    short cur_score {0}, minscore {INF}, minscoreidx {-1};
    
    int s;
    for (s = 0; s <= 2 * (target_len + BAND_SIZE); s += 2) {
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
            cur_score = _mm_extract_epi16(_m1, s / 2 - target_len);
            
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
        
        _d1 = _mm_insert_epi16(_mm_slli_si128(_d1, 2), INF, 0);
        
        _i1 = _mm_add_epi16(_mm_min_epi16(_mm_add_epi16(_i2, _gap_extend),
                                          _mm_add_epi16(_m2, _gap_open)), _nuc_prior);
        
        _backpointers[s] = _mm_or_si128(_mm_or_si128(_mm_and_si128(_three, _m1),
                                                     _mm_slli_epi16(_mm_and_si128(_three, _i1),
                                                                    2 * insert_label)),
                                        _mm_slli_epi16(_mm_and_si128(_three, _d1), 2 * delete_label));
        
        // set state labels
        _m1 = _mm_andnot_si128(_three, _m1);
        _i1 = _mm_or_si128(_mm_andnot_si128(_three, _i1), _mm_srli_epi16(_three, 1));
        _d1 = _mm_or_si128(_mm_andnot_si128(_three, _d1), _three);
        
        // S odd
        
        // truth needs updating; target is current
        const auto pos = BAND_SIZE + s / 2;
        
        const char base {(pos < truth_len) ? truth[pos] : 'N'};
        
        _truthwin   = _mm_insert_epi16(_mm_srli_si128(_truthwin, 2), base, BAND_SIZE - 1);
        
        _truthnqual = _mm_insert_epi16(_mm_srli_si128(_truthnqual, 2),
                                       (base == 'N') ? N_SCORE : INF, BAND_SIZE - 1);
        
        _snvmaskwin   = _mm_insert_epi16(_mm_srli_si128(_snvmaskwin, 2),
                                         pos < truth_len ? snv_mask[pos] : 'N', BAND_SIZE - 1);
        
        _snv_priorwin = _mm_insert_epi16(_mm_srli_si128(_snv_priorwin, 2),
                                         (pos < truth_len) ? snv_prior[pos] << 2 : INF << 2,
                                         BAND_SIZE - 1);
        
        _gap_open  = _mm_insert_epi16(_mm_srli_si128(_gap_open,  2),
                                      gap_open[pos < truth_len ? pos : truth_len - 1] << 2,
                                      BAND_SIZE - 1);
        
        _initmask  = _mm_slli_si128(_initmask, 2);
        _initmask2 = _mm_slli_si128(_initmask2, 2);
        
        _m2 = _mm_min_epi16(_m2, _mm_min_epi16(_i2, _d2));
        
        // at this point, extract minimum score.  Referred-to position must
        // be y==target_len-1, so that current position has y==target_len; i==0 so d=0 and y=s/2
        if (s / 2 >= target_len) {
            cur_score = _mm_extract_epi16(_m2, s / 2 - target_len);
            
            if (cur_score < minscore) {
                minscore = cur_score;
                minscoreidx = s+1;
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
                               INF, BAND_SIZE - 1);
        
        _backpointers[s + 1] = _mm_or_si128(_mm_or_si128(_mm_and_si128(_three, _m2),
                                                         _mm_slli_epi16(_mm_and_si128(_three, _i2), 2 * insert_label)),
                                            _mm_slli_epi16(_mm_and_si128(_three, _d2), 2 * delete_label));
        
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
    auto state  = ((((short*)(_backpointers + s))[i]) >> (2 * match_label)) & 3;
    
    s -= 2;
    
    // this is 2*y (s even) or 2*y+1 (s odd)
    while (y > 0) {
        const auto new_state = ((((short*)(_backpointers + s))[i]) >> (2 * state)) & 3;
        
        if (state == match_label) {
            s -= 2;
            aln1[alnidx] = truth[--x];
            aln2[alnidx] = target[--y];
        } else if (state == insert_label) {
            i += s & 1;
            s -= 1;
            aln1[alnidx] = Gap;
            aln2[alnidx] = target[--y];
        } else {
            s -= 1;
            i -= s & 1;
            aln1[alnidx] = truth[--x];
            aln2[alnidx] = Gap;
        }
        
        state = new_state;
        alnidx++;
    }
    
    aln1[alnidx] = 0;
    aln2[alnidx] = 0;
    
    first_pos = x;
    
    // reverse them
    int j;
    for (i = 0, j = alnidx - 1; i < j; ++i, j--) {
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
                          const int first_pos, const char* aln1, const char* aln2)
{
    static constexpr char MATCH {'M'}, INSERTION {'I'}, DELETION {'D'};
    
    auto prev_state = MATCH;
    
    int x {first_pos}; // index into haplotype
    int y {0};         // index into read
    int i {0};         // index into alignment
    
    int result {0}; // alignment score (within flank)
    
    while (aln1[i]) {
        auto new_state = MATCH;
        
        if (aln1[i] == Gap) new_state = INSERTION;
        if (aln2[i] == Gap) new_state = DELETION;  // can't be both '-'
        
        switch (new_state) {
            case MATCH:
            {
                if ((aln1[i] != aln2[i]) && (x < lhs_flank_len || x >= (truth_len - rhs_flank_len))) {
                    if (aln1[i] != 'N') {
                        result += quals[y];
                    } else {
                        result += N_SCORE >> 2;
                    }
                }
                ++x;
                ++y;
                break;
            }
            case INSERTION:
            {
                if (x < lhs_flank_len || x >= (truth_len - rhs_flank_len)) {
                    if (prev_state == INSERTION) {
                        result += gap_extend + nuc_prior;
                    } else {
                        // gap open score is charged for insertions just after the corresponding base,
                        // hence the -1
                        result += gap_open[x - 1] + nuc_prior;
                    }
                }
                ++y;
                break;
            }
            case DELETION:
            {
                if (x < lhs_flank_len || x >= (truth_len - rhs_flank_len)) {
                    if (prev_state == DELETION) {
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
                          const int first_pos, const char* aln1, const char* aln2)
{
    static constexpr char MATCH {'M'}, INSERTION {'I'}, DELETION {'D'};
    
    auto prev_state = MATCH;
    
    int x {first_pos}; // index into truth
    int y {0};         // index into target
    int i {0};         // index into alignment
    
    int result {0}; // alignment score (within flank)
    
    const auto rhs_flank_begin = truth_len - rhs_flank_len;
    
    while (aln1[i]) {
        auto new_state = MATCH;
        
        if (aln1[i] == Gap) new_state = INSERTION;
        if (aln2[i] == Gap) new_state = DELETION;  // can't be both '-'
        
        switch (new_state) {
            case MATCH:
            {
                if ((aln1[i] != aln2[i]) && (x < lhs_flank_len || x >= rhs_flank_begin)) {
                    if (aln1[i] != 'N') {
                        result += (snv_mask[x] == target[y]) ? std::min(quals[y], snv_prior[x]) : quals[y];
                    } else {
                        result += N_SCORE >> 2;
                    }
                }
                ++x;
                ++y;
                break;
            }
            case INSERTION:
            {
                if (x < lhs_flank_len || x >= rhs_flank_begin) {
                    if (prev_state == INSERTION) {
                        result += gap_extend + nuc_prior;
                    } else {
                        // gap open score is charged for insertions just after the corresponding base,
                        // hence the -1
                        result += gap_open[x - 1] + nuc_prior;
                    }
                }
                ++y;
                break;
            }
            case DELETION:
            {
                if (x < lhs_flank_len || x >= rhs_flank_begin) {
                    if (prev_state == DELETION) {
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
} // namespace SimdPairHmm
