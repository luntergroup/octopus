// Copyright (c) 2015-2021 Daniel Cooke and Gerton Lunter
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef simd_pair_hmm_factory_hpp
#define simd_pair_hmm_factory_hpp

#include <type_traits>

#include "simd_pair_hmm.hpp"
#include "sse2_pair_hmm_impl.hpp"
#include "avx2_pair_hmm_impl.hpp"
#include "avx512_pair_hmm_impl.hpp"
#include "neon_pair_hmm_impl.hpp"
#include "rolling_initializer.hpp"

namespace octopus { namespace hmm { namespace simd {

#if defined(SSE2_PHMM)
template <unsigned BandSize,
          typename ScoreType = short,
          template <class> class InitializerType = InsertRollingInitializer>
using SSE2PairHMM = PairHMM<SSE2PairHMMInstructionSet<BandSize, ScoreType>, InitializerType>;

using FastestSSE2PairHMM = SSE2PairHMM<8, short, InsertRollingInitializer>;

#endif // defined(SSE2_PHMM)

#if defined(AVX512_PHMM)

template <unsigned BandSize,
          typename ScoreType = short,
          template <class> class InitializerType = InsertRollingInitializer>
using AVX512PairHMM = PairHMM<AVX512PairHMMInstructionSet<BandSize, ScoreType>, InitializerType>;

#endif // defined(AVX512_PHMM)

#if defined(AVX2_PHMM)

template <unsigned BandSize,
          typename ScoreType = short,
          template <class> class InitializerType = InsertRollingInitializer>
using AVX2PairHMM = PairHMM<AVX2PairHMMInstructionSet<BandSize, ScoreType>, InitializerType>;

#endif // defined(AVX2_PHMM)

#if defined(NEON_PHMM)

template <unsigned BandSize,
          typename ScoreType = short,
          template <class> class InitializerType = InsertRollingInitializer>
using NEONPairHMM = PairHMM<NEONPairHMMInstructionSet<BandSize, ScoreType>, InitializerType>;

#endif // defined(NEON_PHMM)

namespace detail {

#if defined(AVX512_PHMM)

template <unsigned BandSize, typename ScoreType>
constexpr bool is_viable_avx2 = BandSize % (32 / sizeof(ScoreType)) == 0;

template <unsigned BandSize, typename ScoreType>
constexpr bool is_viable_avx512 = BandSize % (64 / sizeof(ScoreType)) == 0;

template <unsigned BandSize,
          typename ScoreType,
          template <class> class InitializerType>
struct PairHMMSelector :
    public std::conditional<
        is_viable_avx512<BandSize, ScoreType>,
        AVX512PairHMM<BandSize, ScoreType, InitializerType>,
        std::conditional_t<
            is_viable_avx2<BandSize, ScoreType>,
            AVX2PairHMM<BandSize, ScoreType, InitializerType>,
            SSE2PairHMM<BandSize, ScoreType, InitializerType>>
    > {};

#elif defined(AVX2_PHMM)

template <unsigned BandSize, typename ScoreType>
constexpr bool is_viable_avx2 = BandSize % (32 / sizeof(ScoreType)) == 0;

template <unsigned BandSize,
          typename ScoreType,
          template <class> class InitializerType>
struct PairHMMSelector :
    public std::conditional<
        is_viable_avx2<BandSize, ScoreType>,
        AVX2PairHMM<BandSize, ScoreType, InitializerType>,
        SSE2PairHMM<BandSize, ScoreType, InitializerType>
    > {};

#elif defined(SSE2_PHMM)

template <unsigned MinBandSize,
typename ScoreType,
template <class> class InitializerType>
struct PairHMMSelector
{
    using type = SSE2PairHMM<MinBandSize, ScoreType, InitializerType>;
};

#elif defined(NEON_PHMM)

template <unsigned MinBandSize,
typename ScoreType,
template <class> class InitializerType>
struct PairHMMSelector
{
    using type = NEONPairHMM<MinBandSize, ScoreType, InitializerType>;
};

#else

template <unsigned MinBandSize,
typename ScoreType,
template <class> class InitializerType>
struct PairHMMSelector
{
    using type = void;
};

#endif

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
