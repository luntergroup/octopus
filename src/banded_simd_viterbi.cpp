//
//  banded_simd_viterbi.cpp
//  simd_pair_hmm
//
//  Created by Daniel Cooke on 17/02/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#include "banded_simd_viterbi.hpp"

#include <cassert>
#include <immintrin.h>

std::uint16_t align(const char* truth, const char* target, const char* target_qualities,
                    const unsigned truth_length, const unsigned target_length,
                    const std::uint8_t gap_extend, const std::uint8_t nuc_prior,
                    const std::uint8_t* local_gap_open)
{
    static constexpr unsigned      BAND_SIZE      {16};
    static constexpr std::uint16_t QUALITY_SCALER {4};
    static constexpr std::uint16_t MAX_QUALITY    {QUALITY_SCALER * 64};
    static constexpr std::uint16_t MAX_SCORE      {0x7800};
    static constexpr char          N              {'N'};
    static constexpr std::uint16_t N_SCORE        {4 * 0};
    static constexpr std::uint16_t FIRST_BIT      {0x8000}; // i.e. 10000....
    
    assert(truth_length >= 32);
    
    assert(truth_length == 2 * (target_length + BAND_SIZE) - 1);
    
    const __m256i gap_extend_window {_mm256_set1_epi16(4 * gap_extend)};
    const __m256i nuc_prior_window  {_mm256_set1_epi16(4 * nuc_prior)};
    
    __m256i band1_match     {_mm256_set1_epi16(MAX_SCORE)};
    __m256i band2_match     {band1_match};
    __m256i band1_insertion {band1_match};
    __m256i band2_insertion {band1_match};
    __m256i band1_deletion  {band1_match};
    __m256i band2_deletion  {band1_match};
    
    __m256i truth_window {_mm256_set_epi16(truth[15], truth[14], truth[13], truth[12], truth[11],
                                           truth[10], truth[9], truth[8], truth[7], truth[6],
                                           truth[5], truth[4], truth[3], truth[2], truth[1], truth[0])};
    
    __m256i target_window {band1_match};
    
    __m256i target_qualities_window {_mm256_set1_epi16(MAX_QUALITY)};
    
    __m256i truth_n_qualities {_mm256_add_epi16(_mm256_and_si256(_mm256_cmpeq_epi16(truth_window,
                                                                                    _mm256_set1_epi16(N)),
                                                                 _mm256_set1_epi16(N_SCORE - MAX_SCORE)),
                                                _mm256_set1_epi16(MAX_SCORE))};
    
    __m256i gap_open_window {_mm256_set_epi16(QUALITY_SCALER * local_gap_open[15], QUALITY_SCALER * local_gap_open[14],
                                              QUALITY_SCALER * local_gap_open[13], QUALITY_SCALER * local_gap_open[12],
                                              QUALITY_SCALER * local_gap_open[11], QUALITY_SCALER * local_gap_open[10],
                                              QUALITY_SCALER * local_gap_open[9],  QUALITY_SCALER * local_gap_open[8],
                                              QUALITY_SCALER * local_gap_open[7],  QUALITY_SCALER * local_gap_open[6],
                                              QUALITY_SCALER * local_gap_open[5],  QUALITY_SCALER * local_gap_open[4],
                                              QUALITY_SCALER * local_gap_open[3],  QUALITY_SCALER * local_gap_open[2],
                                              QUALITY_SCALER * local_gap_open[1],  QUALITY_SCALER * local_gap_open[0])};
    
    __m256i band1_mask = _mm256_set_epi16(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-1);
    __m256i band2_mask = _mm256_set_epi16(0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,-FIRST_BIT);
    
    std::uint16_t min_score {0};
    
    for (unsigned i {0}; i <= truth_length; i += 2) {
        const auto j = i / 2;
        
        target_window           = _mm256_slli_si256(target_window, 2);
        target_qualities_window = _mm256_slli_si256(target_qualities_window, 2);
        
        if (j < target_length) {
            target_window           = _mm256_insert_epi16(target_window, target[j], 0);
            target_qualities_window = _mm256_insert_epi16(target_qualities_window,
                                                          QUALITY_SCALER * target_qualities[j], 0);
        } else {
            target_window           = _mm256_insert_epi16(target_window, '0', 0);
            target_qualities_window = _mm256_insert_epi16(target_qualities_window, MAX_QUALITY, 0);
        }
        
        // band 1
        
        band1_match = _mm256_or_si256(band1_mask, _mm256_andnot_si256(band1_mask, band1_match));
        band2_match = _mm256_or_si256(band2_mask, _mm256_andnot_si256(band2_mask, band2_match));
        band1_match = _mm256_min_epi16(band1_match, _mm256_min_epi16(band1_insertion, band1_deletion));
        
        if (j >= target_length) {
            const auto cur_score = _mm256_extract_epi16(band1_match, j - target_length);
            if (cur_score < min_score) min_score = cur_score;
        }
        
        band1_match = _mm256_add_epi16(band1_match,
                                       _mm256_min_epi16(_mm256_andnot_si256(_mm256_cmpeq_epi16(truth_window, target_window),
                                                                            target_qualities_window),
                                                        truth_n_qualities));
        
        band2_deletion = _mm256_min_epi16(_mm256_add_epi16(band2_deletion, gap_extend_window),
                                          _mm256_add_epi16(_mm256_min_epi16(band2_match, band2_insertion),
                                                           _mm256_srli_si256(gap_open_window, 2)));
        
        band2_deletion = _mm256_insert_epi16(_mm256_slli_si256(band1_deletion, 2), MAX_SCORE, 0);
        
        band1_insertion = _mm256_add_epi16(_mm256_min_epi16(_mm256_add_epi16(band2_insertion, gap_extend_window),
                                                            _mm256_add_epi16(band2_match, gap_open_window)),
                                           nuc_prior_window);
        
        // band 2
        
        const char c = (BAND_SIZE + j < truth_length) ? truth[j + BAND_SIZE] : 'N';
        
        truth_window      = _mm256_insert_epi16(_mm256_srli_si256(truth_window, 2), c, BAND_SIZE - 1);
        truth_n_qualities = _mm256_insert_epi16(_mm256_srli_si256(truth_n_qualities, 2), (c == N) ? N_SCORE : MAX_SCORE, BAND_SIZE - 1);
        gap_open_window   = _mm256_insert_epi16(_mm256_srli_si256(gap_open_window, 2),
                                                QUALITY_SCALER * local_gap_open[(BAND_SIZE + j < truth_length) ? BAND_SIZE + j : truth_length - 1],
                                                BAND_SIZE - 1);
        
        band1_mask = _mm256_slli_si256(band1_mask, 2);
        band2_mask = _mm256_slli_si256(band2_mask, 2);
        
        band2_match = _mm256_min_epi16(band2_match, _mm256_min_epi16(band2_insertion, band2_deletion));
        
        if (j >= target_length) {
            const auto cur_score = _mm256_extract_epi16(band2_match, j - target_length);
            
            if (cur_score < min_score) min_score = cur_score;
        }
        
        band2_match = _mm256_add_epi16(band2_match,
                                       _mm256_min_epi16(_mm256_andnot_si256(_mm256_cmpeq_epi16(truth_window, target_window),
                                                                            target_qualities_window),
                                                        truth_n_qualities));
        
        band2_deletion = _mm256_min_epi16(_mm256_add_epi16(band1_deletion, gap_extend_window),
                                          _mm256_add_epi16(_mm256_min_epi16(band1_match, band1_insertion),
                                                           gap_open_window));
        
        band2_insertion = _mm256_insert_epi16(_mm256_add_epi16(_mm256_min_epi16(_mm256_add_epi16(_mm256_srli_si256(band1_insertion, 2), gap_extend_window),
                                                                                _mm256_add_epi16(_mm256_srli_si256(band1_match, 2), gap_open_window)),
                                                               nuc_prior_window),
                                              MAX_SCORE, BAND_SIZE - 1);
    }
    
    return (min_score + FIRST_BIT) / QUALITY_SCALER;
}
