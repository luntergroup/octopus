// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef hiseq_indel_error_model_hpp
#define hiseq_indel_error_model_hpp

#include "indel_error_model.hpp"
#include "repeat_based_indel_error_model.hpp"

namespace octopus {

class HiSeqIndelErrorModel : public RepeatBasedIndelErrorModel
{
public:
    using IndelErrorModel::PenaltyType;
    using IndelErrorModel::PenaltyVector;
    
    HiSeqIndelErrorModel() = default;
    
    HiSeqIndelErrorModel(const HiSeqIndelErrorModel&)            = default;
    HiSeqIndelErrorModel& operator=(const HiSeqIndelErrorModel&) = default;
    HiSeqIndelErrorModel(HiSeqIndelErrorModel&&)                 = default;
    HiSeqIndelErrorModel& operator=(HiSeqIndelErrorModel&&)      = default;
    
    ~HiSeqIndelErrorModel() = default;
    
private:
    virtual std::unique_ptr<IndelErrorModel> do_clone() const override;
    virtual PenaltyType get_default_open_penalty() const noexcept override;
    virtual PenaltyType get_open_penalty(const Sequence& motif, unsigned length) const noexcept override;
    virtual PenaltyType get_default_extension_penalty() const noexcept override;
    virtual PenaltyType get_extension_penalty(const Sequence& motif, unsigned length) const noexcept override;
};

} // namespace octopus

#endif
