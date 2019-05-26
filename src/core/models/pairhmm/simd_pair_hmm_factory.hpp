// Copyright (c) 2015-2019 Daniel Cooke and Gerton Lunter
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef simd_pair_hmm_factory_hpp
#define simd_pair_hmm_factory_hpp

#include <type_traits>
#include <tuple>

#include <boost/variant.hpp>

#include "simd_pair_hmm.hpp"
#include "sse2_pair_hmm_impl.hpp"
#ifdef __AVX2__
    #include "avx2_pair_hmm_impl.hpp"
#endif
#ifdef __AVX512__
    #include "avx512_pair_hmm_impl.hpp"
#endif /* __AVX2__ */

namespace octopus { namespace hmm { namespace simd {

template <unsigned MinBandSize,
          typename ScoreType = short,
          typename RollingInitializer = InsertRollingInitializer<SSE2PairHMMInstructionSet<MinBandSize, ScoreType>>>
using SSE2PairHMM = PairHMM<SSE2PairHMMInstructionSet<MinBandSize, ScoreType>, RollingInitializer>;

using FastestSSE2PairHMM = SSE2PairHMM<8, short>;

#ifdef __AVX2__

using AVX2PairHMM = PairHMM<AVX2PairHMMInstructionSet>;

#endif /* __AVX2__ */

namespace detail {

#ifdef __AVX2__

template <unsigned MinBandSize, typename ScoreType>
struct PairHMMSelector:
    public std::conditional<
        MinBandSize < AVX2PairHMM::band_size(),
        SSE2PairHMM<MinBandSize, ScoreType>,
        AVX2PairHMM
    > {};

#else

template <unsigned MinBandSize, typename ScoreType>
struct PairHMMSelector
{
    using type = SSE2PairHMM<MinBandSize, ScoreType>;
};

#endif /* __AVX2__ */

} // namespace detail

template <unsigned MinBandSize, typename ScoreType = short>
using SimdPairHMM = typename detail::PairHMMSelector<MinBandSize, ScoreType>::type;

template <unsigned MinBandSize, typename ScoreType = short>
auto make_simd_pair_hmm() noexcept { return SimdPairHMM<MinBandSize, ScoreType> {}; }

namespace detail {

using PairHMMs = std::tuple<
    SimdPairHMM<8>,
    SimdPairHMM<16>,
    SimdPairHMM<24>,
    SimdPairHMM<32>
    >;

using PairHmmVariants = boost::variant<
    std::tuple_element_t<0, PairHMMs>,
    std::tuple_element_t<1, PairHMMs>,
    std::tuple_element_t<2, PairHMMs>
    >;

const static PairHMMs pair_hmms = {};

} // namespace detail

detail::PairHmmVariants make_simd_pair_hmm(const unsigned band_size) noexcept
{
    if (band_size <= std::tuple_element_t<0, detail::PairHMMs>::band_size()) {
        return std::get<0>(detail::pair_hmms);
    } else if (band_size <= std::tuple_element_t<1, detail::PairHMMs>::band_size()) {
        return std::get<1>(detail::pair_hmms);
    } else {
        return std::get<2>(detail::pair_hmms);
    }
}

} // namespace simd
} // namespace hmm
} // namespace octopus

#endif
