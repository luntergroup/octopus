// Copyright (c) 2015-2019 Daniel Cooke and Gerton Lunter
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef simd_pair_hmm_wrapper_hpp
#define simd_pair_hmm_wrapper_hpp

#include <tuple>

#include <boost/variant.hpp>

namespace octopus { namespace hmm { namespace simd {

class PairHMMWrapper
{
public:
    PairHMMWrapper(unsigned band_size = 8) noexcept { reset(band_size); }
    
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
    
    void reset(unsigned band_size) noexcept
    {
        hmm_ = make_simd_pair_hmm(band_size);
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
    using PairHMMs = std::tuple<
        SimdPairHMM<8,  short>,
        SimdPairHMM<16, short>,
        SimdPairHMM<24, short>,
        SimdPairHMM<32, short>,
        SimdPairHMM<40, short>,
        SimdPairHMM<48, short>,
        SimdPairHMM<4,  int>,
        SimdPairHMM<8,  int>,
        SimdPairHMM<12, int>,
        SimdPairHMM<16, int>,
        SimdPairHMM<20, int>,
        SimdPairHMM<24, int>,
        SimdPairHMM<28, int>,
        SimdPairHMM<32, int>
    >;
    
    using PairHmmVariants = boost::variant<
        std::tuple_element_t<0, PairHMMs>,
        std::tuple_element_t<1, PairHMMs>,
        std::tuple_element_t<2, PairHMMs>
    >;
    
    PairHmmVariants hmm_;
    
    PairHmmVariants make_simd_pair_hmm(const unsigned band_size) noexcept
    {
        const static PairHMMs pair_hmms {};
        if (band_size <= std::tuple_element_t<0, PairHMMs>::band_size()) {
            return std::get<0>(pair_hmms);
        } else if (band_size <= std::tuple_element_t<1, PairHMMs>::band_size()) {
            return std::get<1>(pair_hmms);
        } else {
            return std::get<2>(pair_hmms);
        }
    }
};

} // namespace simd
} // namespace hmm
} // namespace octopus

#endif
