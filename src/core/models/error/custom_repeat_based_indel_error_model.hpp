// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef custom_repeat_based_indel_error_model_hpp
#define custom_repeat_based_indel_error_model_hpp

#include <unordered_map>
#include <string>
#include <vector>

#include <boost/optional.hpp>

#include "indel_error_model.hpp"
#include "repeat_based_indel_error_model.hpp"

namespace octopus {

class CustomRepeatBasedIndelErrorModel : public RepeatBasedIndelErrorModel
{
public:
    using IndelErrorModel::PenaltyType;
    using IndelErrorModel::PenaltyVector;
    
    using MotifPenaltyMap = std::unordered_map<Sequence, PenaltyVector>;
    
    CustomRepeatBasedIndelErrorModel() = delete;
    
    CustomRepeatBasedIndelErrorModel(MotifPenaltyMap gap_open_penalties, PenaltyType extend_penalty);
    CustomRepeatBasedIndelErrorModel(MotifPenaltyMap gap_open_penalties, MotifPenaltyMap gap_extend_penalties);
    
    CustomRepeatBasedIndelErrorModel(const CustomRepeatBasedIndelErrorModel&)            = default;
    CustomRepeatBasedIndelErrorModel& operator=(const CustomRepeatBasedIndelErrorModel&) = default;
    CustomRepeatBasedIndelErrorModel(CustomRepeatBasedIndelErrorModel&&)                 = default;
    CustomRepeatBasedIndelErrorModel& operator=(CustomRepeatBasedIndelErrorModel&&)      = default;
    
    ~CustomRepeatBasedIndelErrorModel() = default;

private:
    MotifPenaltyMap gap_open_penalties_;
    boost::optional<MotifPenaltyMap> gap_extend_penalties_;
    PenaltyType default_gap_open_, default_gap_extend_;
    std::vector<Sequence> ns_;
    
    std::unique_ptr<IndelErrorModel> do_clone() const override;
    PenaltyType get_default_open_penalty() const noexcept override;
    PenaltyType get_open_penalty(const Sequence& motif, unsigned length) const noexcept override;
    PenaltyType get_default_extension_penalty() const noexcept override;
    PenaltyType get_extension_penalty(const Sequence& motif, unsigned length) const noexcept override;
};

CustomRepeatBasedIndelErrorModel::MotifPenaltyMap make_penalty_map(std::string model);

} // namespace octopus

#endif
