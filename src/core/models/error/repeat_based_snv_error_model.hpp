// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef repeat_based_snv_error_model_hpp
#define repeat_based_snv_error_model_hpp

#include <vector>
#include <array>
#include <cstdint>

#include "snv_error_model.hpp"

namespace octopus {

class Haplotype;

class BasicRepeatBasedSNVErrorModel : public SnvErrorModel
{
public:
    using SnvErrorModel::MutationVector;
    using SnvErrorModel::PenaltyType;
    using SnvErrorModel::PenaltyVector;
    
    struct Parameters
    {
        PenaltyVector homopolymer_penalty_caps;
        PenaltyVector dinucleotide_penalty_caps;
        PenaltyVector trinucleotide_penalty_caps;
    };
    
    BasicRepeatBasedSNVErrorModel() = delete;
    
    BasicRepeatBasedSNVErrorModel(Parameters params);
    
    BasicRepeatBasedSNVErrorModel(const BasicRepeatBasedSNVErrorModel&)            = default;
    BasicRepeatBasedSNVErrorModel& operator=(const BasicRepeatBasedSNVErrorModel&) = default;
    BasicRepeatBasedSNVErrorModel(BasicRepeatBasedSNVErrorModel&&)                 = default;
    BasicRepeatBasedSNVErrorModel& operator=(BasicRepeatBasedSNVErrorModel&&)      = default;
    
    virtual ~BasicRepeatBasedSNVErrorModel() = default;

private:
    static constexpr std::size_t max_period_ = 3;
    std::array<std::array<PenaltyType, 51>, max_period_> penalty_caps_;
    
    virtual std::unique_ptr<SnvErrorModel> do_clone() const override;
    virtual void do_evaluate(const Haplotype& haplotype,
                             MutationVector& forward_snv_mask, PenaltyVector& forward_snv_priors,
                             MutationVector& reverse_snv_mask, PenaltyVector& reverse_snv_priors) const override ;
};

} // namespace octopus

#endif
