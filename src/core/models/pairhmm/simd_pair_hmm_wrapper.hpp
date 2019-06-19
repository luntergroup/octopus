// Copyright (c) 2015-2019 Daniel Cooke and Gerton Lunter
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef simd_pair_hmm_wrapper_hpp
#define simd_pair_hmm_wrapper_hpp

#include <tuple>

#include <boost/variant.hpp>

#include "exceptions/user_error.hpp"
#include "simd_pair_hmm_factory.hpp"

namespace octopus { namespace hmm { namespace simd {

namespace detail {

template <typename T, unsigned step,  std::size_t... I>
constexpr auto make_phmm_tuple(std::index_sequence<I...>) { return std::make_tuple(SimdPairHMM<step * (I + 1), T>{}...); }

template <typename Tuple>
struct to_variant;

template <typename... Ts>
struct to_variant<std::tuple<Ts...>>
{
    using type = boost::variant<Ts...>;
};

template <typename Tuple>
using to_variant_t = typename to_variant<Tuple>::type;

} // namespace detail

class PairHMMWrapper
{
public:
    enum class ScorePrecision { int16, int32 };
    
    class TooLargeBandSizeError : public std::runtime_error
    {
    public:
        TooLargeBandSizeError() = delete;
        
        TooLargeBandSizeError(int requested, int max)
        : std::runtime_error {"requested band size is too large"}
        , requested_ {requested}
        , max_ {max}
        {}
    
        int requested() const noexcept { return requested_; };
        int max() const noexcept { return max_; }
    
    private:
        int requested_, max_;
    };
    
    PairHMMWrapper(int min_band_size = 8, ScorePrecision score_precision = ScorePrecision::int16)
    {
        reset(min_band_size, score_precision);
    }
    
    PairHMMWrapper(const PairHMMWrapper&)            = default;
    PairHMMWrapper& operator=(const PairHMMWrapper&) = default;
    PairHMMWrapper(PairHMMWrapper&&)                 = default;
    PairHMMWrapper& operator=(PairHMMWrapper&&)      = default;
    
    ~PairHMMWrapper() = default;
    
    int band_size() const noexcept
    {
        return boost::apply_visitor([] (const auto& hmm) noexcept { return hmm.band_size(); }, hmm_);
    }
    
    const char* name() const noexcept
    {
        return boost::apply_visitor([] (const auto& hmm) noexcept { return hmm.name(); }, hmm_);
    }
    
    void reset(int min_band_size, ScorePrecision score_precision = ScorePrecision::int16)
    {
        hmm_ = make_simd_pair_hmm(min_band_size, score_precision);
    }
    
    template <typename OpenPenaltyArrayOrConstant,
              typename ExtendPenaltyArrayOrConstant>
    int
    align(const char* truth,
          const char* target,
          const std::int8_t* qualities,
          const int truth_len,
          const int target_len,
          const OpenPenaltyArrayOrConstant gap_open,
          const ExtendPenaltyArrayOrConstant gap_extend,
          short nuc_prior) const noexcept
    {
        return boost::apply_visitor([&] (const auto& hmm) noexcept {
            return hmm.align(truth, target, qualities, truth_len, target_len, gap_open, gap_extend, nuc_prior);
        }, hmm_);
    }
    template <typename OpenPenaltyArrayOrConstant,
              typename ExtendPenaltyArrayOrConstant>
    int
    align(const char* truth,
          const char* target,
          const std::int8_t* qualities,
          const int truth_len,
          const int target_len,
          const char* snv_mask,
          const std::int8_t* snv_prior,
          const OpenPenaltyArrayOrConstant gap_open,
          const ExtendPenaltyArrayOrConstant gap_extend,
          short nuc_prior) const noexcept
    {
        return boost::apply_visitor([&] (const auto& hmm) noexcept {
            return hmm.align(truth, target, qualities, truth_len, target_len, snv_mask, snv_prior, gap_open, gap_extend, nuc_prior);
        }, hmm_);
    }
    template <typename OpenPenaltyArrayOrConstant,
              typename ExtendPenaltyArrayOrConstant>
    int
    align(const char* truth,
          const char* target,
          const std::int8_t* qualities,
          const int truth_len,
          const int target_len,
          const OpenPenaltyArrayOrConstant gap_open,
          const ExtendPenaltyArrayOrConstant gap_extend,
          short nuc_prior,
          int& first_pos,
          char* align1,
          char* align2) const noexcept
    {
        return boost::apply_visitor([&] (const auto& hmm) noexcept {
            return hmm.align(truth, target, qualities, truth_len, target_len, gap_open, gap_extend, nuc_prior, first_pos, align1, align2);
        }, hmm_);
    }
    template <typename OpenPenaltyArrayOrConstant,
              typename ExtendPenaltyArrayOrConstant>
    int
    align(const char* truth,
          const char* target,
          const std::int8_t* qualities,
          const int truth_len,
          const int target_len,
          const char* snv_mask,
          const std::int8_t* snv_prior,
          const OpenPenaltyArrayOrConstant gap_open,
          const ExtendPenaltyArrayOrConstant gap_extend,
          short nuc_prior,
          int& first_pos,
          char* align1,
          char* align2) const noexcept
    {
        return boost::apply_visitor([&] (const auto& hmm) noexcept {
            return hmm.align(truth, target, qualities, truth_len, target_len, snv_mask, snv_prior, gap_open, gap_extend, nuc_prior, first_pos, align1, align2);
        }, hmm_);
    }
    template <typename OpenPenaltyArrayOrConstant,
              typename ExtendPenaltyArrayOrConstant>
    int
    calculate_flank_score(const int truth_len,
                          const int lhs_flank_len,
                          const int rhs_flank_len,
                          const std::int8_t* quals,
                          const OpenPenaltyArrayOrConstant gap_open,
                          const ExtendPenaltyArrayOrConstant gap_extend,
                          const short nuc_prior,
                          const int first_pos,
                          const char* aln1,
                          const char* aln2,
                          int& target_mask_size) const noexcept
    {
        return boost::apply_visitor([&] (const auto& hmm) noexcept {
            return hmm.calculate_flank_score(truth_len, lhs_flank_len, rhs_flank_len, quals, gap_open, gap_extend, nuc_prior, first_pos, aln1, aln2, target_mask_size);
        }, hmm_);
    }
    template <typename OpenPenaltyArrayOrConstant,
              typename ExtendPenaltyArrayOrConstant>
    int
    calculate_flank_score(int truth_len,
                          int lhs_flank_len,
                          int rhs_flank_len,
                          const char* target,
                          const std::int8_t* quals,
                          const char* snv_mask,
                          const std::int8_t* snv_prior,
                          const OpenPenaltyArrayOrConstant gap_open,
                          const ExtendPenaltyArrayOrConstant gap_extend,
                          short nuc_prior,
                          int first_pos,
                          const char* aln1,
                          const char* aln2,
                          int& target_mask_size) const noexcept
    {
        return boost::apply_visitor([&] (const auto& hmm) noexcept {
            return hmm.calculate_flank_score(truth_len, lhs_flank_len, rhs_flank_len, target, quals, snv_mask, snv_prior, gap_open, gap_extend, nuc_prior, first_pos, aln1, aln2, target_mask_size);
        }, hmm_);
    }

private:
    using ShortPairHMMs = decltype(detail::make_phmm_tuple<short, 8>(std::make_index_sequence<10>()));
    using IntPairHMMs   = decltype(detail::make_phmm_tuple<int, 8>(std::make_index_sequence<10>()));
    using PairHMMs      = decltype(std::tuple_cat(ShortPairHMMs {}, IntPairHMMs {}));
    
    using PairHmmVariant = detail::to_variant_t<PairHMMs>;
    
    PairHmmVariant hmm_;
    
    template <typename Hmms>
    constexpr static int max_band_size() { return std::tuple_element_t<std::tuple_size<Hmms>::value - 1, Hmms>::band_size(); }
    
    template <typename Hmms, std::size_t... Is>
    PairHmmVariant make_simd_pair_hmm_helper(const int min_band_size, std::index_sequence<Is...>) const
    {
        if (min_band_size <= max_band_size<Hmms>()) {
            PairHmmVariant result {};
            int unused[] = {((Is == 0 || std::tuple_element_t<(Is > 0 ? Is - 1 : 0), Hmms>::band_size() < min_band_size)
                            && min_band_size <= std::tuple_element_t<Is, Hmms>::band_size()
                            ? (result = std::tuple_element_t<Is, Hmms> {}, 0) : 0)...};
            (void) unused;
            return result;
        } else {
            throw TooLargeBandSizeError {min_band_size, max_band_size<Hmms>()};
        }
    }
    PairHmmVariant make_simd_pair_hmm(const int min_band_size, const ScorePrecision score_precision) const
    {
        if (score_precision == ScorePrecision::int16) {
            using Hmms = ShortPairHMMs;
            return make_simd_pair_hmm_helper<Hmms>(min_band_size, std::make_index_sequence<std::tuple_size<Hmms>::value>());
        } else {
            using Hmms = IntPairHMMs;
            return make_simd_pair_hmm_helper<Hmms>(min_band_size, std::make_index_sequence<std::tuple_size<Hmms>::value>());
        }
    }

public:
    static int max_band_size(ScorePrecision score_precision) noexcept
    {
        if (score_precision == ScorePrecision::int16) {
            return max_band_size<ShortPairHMMs>();
        } else {
            return max_band_size<IntPairHMMs>();
        }
    }
};

} // namespace simd
} // namespace hmm
} // namespace octopus

#endif
