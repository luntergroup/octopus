// Copyright (c) 2015-2019 Daniel Cooke and Gerton Lunter
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef simd_pair_hmm_factory_hpp
#define simd_pair_hmm_factory_hpp

#include <type_traits>

#include "simd_pair_hmm.hpp"
#include "sse2_pair_hmm_impl.hpp"
#include "avx2_pair_hmm_impl.hpp"
//#include "avx512_pair_hmm_impl.hpp"
#include "rolling_initializer.hpp"

namespace octopus { namespace hmm { namespace simd {

template <unsigned BandSize,
          typename ScoreType = short,
          template <class> class InitializerType = InsertRollingInitializer>
using SSE2PairHMM = PairHMM<SSE2PairHMMInstructionSet<BandSize, ScoreType>, InitializerType>;

using FastestSSE2PairHMM = SSE2PairHMM<8, short, InsertRollingInitializer>;

#if defined(AVX2_PHMM)
template <template <class> class InitializerType = InsertRollingInitializer>
using AVX2PairHMM = PairHMM<AVX2PairHMMInstructionSet, InitializerType>;
#endif

namespace detail {

#if defined(AVX2_PHMM)

template <unsigned BandSize,
          typename ScoreType,
          template <class> class InitializerType>
struct PairHMMSelector :
    public std::conditional<
        BandSize == AVX2PairHMM<InitializerType>::band_size() && std::is_same<ScoreType, short>::value,
        AVX2PairHMM<InitializerType>,
        SSE2PairHMM<BandSize, ScoreType, InitializerType>
    > {};

#else // defined(AVX2_PHMM)

template <unsigned MinBandSize,
typename ScoreType,
template <class> class InitializerType>
struct PairHMMSelector
{
    using type = SSE2PairHMM<MinBandSize, ScoreType, InitializerType>;
};

#endif // defined(AVX2_PHMM)

} // namespace detail

template <unsigned BandSize,
          typename ScoreType = short,
          template <class> class InitializerType = InsertRollingInitializer>
using SimdPairHMM = typename detail::PairHMMSelector<BandSize, ScoreType, InitializerType>::type;

template <unsigned BandSize,
          typename ScoreType = short,
          template <class> class InitializerType = InsertRollingInitializer>
constexpr auto make_simd_pair_hmm() { return SimdPairHMM<BandSize, ScoreType, InitializerType> {}; }

constexpr auto make_fastest_simd_pair_hmm() { return make_simd_pair_hmm<8, short, InsertRollingInitializer>(); }

} // namespace simd
} // namespace hmm
} // namespace octopus

#endif
