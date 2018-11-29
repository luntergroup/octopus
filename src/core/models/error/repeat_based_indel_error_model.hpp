// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef repeat_based_indel_error_model_hpp
#define repeat_based_indel_error_model_hpp

#include "indel_error_model.hpp"

#include "core/types/haplotype.hpp"

namespace octopus {

class RepeatBasedIndelErrorModel : public IndelErrorModel
{
public:
    using IndelErrorModel::PenaltyType;
    using IndelErrorModel::PenaltyVector;
    
    RepeatBasedIndelErrorModel() = default;
    
    RepeatBasedIndelErrorModel(const RepeatBasedIndelErrorModel&)            = default;
    RepeatBasedIndelErrorModel& operator=(const RepeatBasedIndelErrorModel&) = default;
    RepeatBasedIndelErrorModel(RepeatBasedIndelErrorModel&&)                 = default;
    RepeatBasedIndelErrorModel& operator=(RepeatBasedIndelErrorModel&&)      = default;

    virtual ~RepeatBasedIndelErrorModel() = default;

protected:
    using Sequence = Haplotype::NucleotideSequence;
    
private:
    void do_set_penalties(const Haplotype& haplotype, PenaltyVector& gap_open_penalties, PenaltyType& gap_extend_penalty) const override;
    void do_set_penalties(const Haplotype& haplotype, PenaltyVector& gap_open_penalties, PenaltyVector& gap_extend_penalties) const override;
    
    virtual std::unique_ptr<IndelErrorModel> do_clone() const override = 0;
    virtual PenaltyType get_default_open_penalty() const noexcept = 0;
    virtual PenaltyType get_open_penalty(const Sequence& motif, unsigned length) const noexcept = 0;
    virtual PenaltyType get_default_extension_penalty() const noexcept = 0;
    virtual PenaltyType get_extension_penalty(const Sequence& motif, unsigned length) const noexcept = 0;
};

} // namespace octopus

#endif
