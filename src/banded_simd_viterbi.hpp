//
//  banded_simd_viterbi.hpp
//  banded_simd_pair_hmm
//
//  Created by Daniel Cooke on 16/02/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#ifndef banded_simd_viterbi_hpp
#define banded_simd_viterbi_hpp

#include <cstdint>

std::uint16_t align(const char* truth, const char* target, const char* target_qualities,
                    unsigned truth_length, unsigned target_length,
                    std::uint8_t gap_extend, std::uint8_t nuc_prior,
                    const std::uint8_t* local_gap_open);

#endif /* banded_simd_viterbi_hpp */
