// Copyright (c) 2015-2019 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef basic_indel_error_model_hpp
#define basic_indel_error_model_hpp

#include <array>

#include "indel_error_model.hpp"
#include "repeat_based_indel_error_model.hpp"

namespace octopus {

class BasicRepeatBasedIndelErrorModel : public RepeatBasedIndelErrorModel
{
public:
    using IndelErrorModel::PenaltyType;
    using IndelErrorModel::PenaltyVector;
    
    struct Parameters
    {
        PenaltyVector AT_homopolymer_open_penalities;
        PenaltyVector CG_homopolymer_open_penalties;
        PenaltyVector dinucleotide_repeat_open_penalties;
        PenaltyVector trinucleotide_repeat_open_penalties;
        PenaltyVector homopolymer_extend_penalties = {3, 3, 3, 3, 3, 3, 4, 5, 6, 6, 8, 8, 7, 6, 5, 4, 3};
        PenaltyVector dinucleotide_repeat_extend_penalties = {3, 3, 5, 4, 3, 2};
        PenaltyVector trinucleotide_repeat_extend_penalties = {3, 3, 5, 4, 3, 2};
    };
    
    BasicRepeatBasedIndelErrorModel() = delete;
    
    BasicRepeatBasedIndelErrorModel(Parameters params);
    
    BasicRepeatBasedIndelErrorModel(const BasicRepeatBasedIndelErrorModel&)            = default;
    BasicRepeatBasedIndelErrorModel& operator=(const BasicRepeatBasedIndelErrorModel&) = default;
    BasicRepeatBasedIndelErrorModel(BasicRepeatBasedIndelErrorModel&&)                 = default;
    BasicRepeatBasedIndelErrorModel& operator=(BasicRepeatBasedIndelErrorModel&&)      = default;
    
    ~BasicRepeatBasedIndelErrorModel() = default;

private:
    using PenaltyArray = std::array<PenaltyType, 50>;
    
    PenaltyArray AT_homopolymer_open_penalities_, CG_homopolymer_open_penalties_,
                    dinucleotide_repeat_open_penalties_, trinucleotide_repeat_open_penalties_;
    PenaltyArray homopolymer_extend_penalties_, dinucleotide_repeat_extend_penalties_, trinucleotide_repeat_extend_penalties_;
    PenaltyType complex_open_penalty_, complex_extend_penalty_;
    
    virtual std::unique_ptr<IndelErrorModel> do_clone() const override;
    virtual PenaltyType get_default_open_penalty() const noexcept override;
    virtual PenaltyType get_open_penalty(const Sequence& motif, unsigned length) const noexcept override;
    virtual PenaltyType get_default_extension_penalty() const noexcept override;
    virtual PenaltyType get_extension_penalty(const Sequence& motif, unsigned length) const noexcept override;
};
    
} // namespace octopus

#endif
