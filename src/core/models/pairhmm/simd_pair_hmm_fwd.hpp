// Copyright (c) 2015-2019 Daniel Cooke and Gerton Lunter
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef simd_pair_hmm_fwd_hpp
#define simd_pair_hmm_fwd_hpp

#include "simd_pair_hmm.hpp"
#include "sse2_pair_hmm_impl.hpp"
#ifdef __AVX2__
    #include "avx2_pair_hmm_impl.hpp"
#endif /* __AVX2__ */

#endif
