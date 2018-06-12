// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef hiseq_indel_error_model_hpp
#define hiseq_indel_error_model_hpp

#include "indel_error_model.hpp"

namespace octopus {

class Haplotype;

class HiSeqIndelErrorModel : public IndelErrorModel
{
public:
    using IndelErrorModel::PenaltyType;
    using IndelErrorModel::PenaltyVector;
    
    HiSeqIndelErrorModel() = default;
    
    HiSeqIndelErrorModel(const HiSeqIndelErrorModel&)            = default;
    HiSeqIndelErrorModel& operator=(const HiSeqIndelErrorModel&) = default;
    HiSeqIndelErrorModel(HiSeqIndelErrorModel&&)                 = default;
    HiSeqIndelErrorModel& operator=(HiSeqIndelErrorModel&&)      = default;
    
private:
    static constexpr std::array<PenaltyType, 50> homopolymerErrors_ =
    {{
     60,60,50,45,41,36,30,25,22,20,19,17,16,15,14,13,12,11,11,10,
     9,9,8,8,7,7,7,6,6,6,5,5,5,4,4,4,3,3,3,3,2,2,2,2,2,1,1,1,1,1
     }};
    static constexpr std::array<PenaltyType, 50> diNucleotideTandemRepeatErrors_ =
    {{
     60,60,48,45,43,41,39,35,31,28,25,21,19,17,15,13,12,11,11,10,
     9,9,8,8,7,7,7,6,6,6,5,5,5,4,4,4,3,3,3,3,2,2,2,2,2,1,1,1,1,1
     }};
    static constexpr std::array<PenaltyType, 50> triNucleotideTandemRepeatErrors_ =
    {{
     60,60,50,48,46,45,42,39,35,31,28,25,22,20,16,14,13,12,12,11,
     10,9,8,8,7,7,7,6,6,6,5,5,5,4,4,4,3,3,3,3,2,2,2,2,2,1,1,1,1,1
     }};
    
    static constexpr std::array<PenaltyType, 50> polyNucleotideTandemRepeatErrors_ =
    {{
     60,60,51,45,45,45,45,45,23,20,19,17,16,15,14,13,12,11,11,10,
     9,9,8,8,7,7,7,6,6,6,5,5,5,4,4,4,3,3,3,3,2,2,2,2,2,1,1,1,1,1
     }};
    static constexpr PenaltyType defaultGapExtension_ = 3;
    
    virtual std::unique_ptr<IndelErrorModel> do_clone() const override;
    virtual PenaltyType do_evaluate(const Haplotype& haplotype, PenaltyVector& gap_open_penalties) const override;
};

} // namespace octopus

#endif
