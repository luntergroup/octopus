// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef bgiseq_indel_error_model_hpp
#define bgiseq_indel_error_model_hpp

#include "indel_error_model.hpp"
#include "repeat_based_indel_error_model.hpp"

namespace octopus {

class BGISeqIndelErrorModel : public RepeatBasedIndelErrorModel
{
public:
    using IndelErrorModel::PenaltyType;
    using IndelErrorModel::PenaltyVector;
    
    BGISeqIndelErrorModel() = default;
    
    BGISeqIndelErrorModel(const BGISeqIndelErrorModel&)            = default;
    BGISeqIndelErrorModel& operator=(const BGISeqIndelErrorModel&) = default;
    BGISeqIndelErrorModel(BGISeqIndelErrorModel&&)                 = default;
    BGISeqIndelErrorModel& operator=(BGISeqIndelErrorModel&&)      = default;
    
    ~BGISeqIndelErrorModel() = default;

private:
    virtual std::unique_ptr<IndelErrorModel> do_clone() const override;
    virtual PenaltyType get_default_open_penalty() const noexcept override;
    virtual PenaltyType get_open_penalty(const Sequence& motif, unsigned length) const noexcept override;
    virtual PenaltyType get_default_extension_penalty() const noexcept override;
    virtual PenaltyType get_extension_penalty(const Sequence& motif, unsigned length) const noexcept override;
};
    
} // namespace octopus

#endif
