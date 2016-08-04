//
//  read_indel_error_model.hpp
//  octopus
//
//  Created by Daniel Cooke on 13/06/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

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
    
    PenaltyType evaluate(const Haplotype& haplotype, PenaltyVector& gap_open_penalties) const;
    
private:
    static constexpr std::array<PenaltyType, 50> Homopolymer_errors_ =
    {{
        45,42,41,39,37,32,28,23,20,19,17,16,15,14,13,12,11,11,10,
        9,9,8,8,7,7,7,6,6,6,5,5,5,4,4,4,3,3,3,3,2,2,2,2,2,1,1,1,1,1,1
    }};
    
    static constexpr std::array<PenaltyType, 50> Di_nucleotide_tandem_repeat_errors_ =
    {{
        45,45,45,45,37,32,28,23,20,19,17,16,15,14,13,12,11,11,10,
        9,9,8,8,7,7,7,6,6,6,5,5,5,4,4,4,3,3,3,3,2,2,2,2,2,1,1,1,1,1,1
    }};
    
    static constexpr std::array<PenaltyType, 50> Tri_nucleotide_tandem_repeat_errors_ =
    {{
        45,45,45,45,45,45,28,23,20,19,17,16,15,14,13,12,11,11,10,
        9,9,8,8,7,7,7,6,6,6,5,5,5,4,4,4,3,3,3,3,2,2,2,2,2,1,1,1,1,1,1
    }};
    
    static constexpr std::array<PenaltyType, 50> Poly_nucleotide_tandem_repeat_errors_ =
    {{
        45,45,45,45,45,45,45,23,20,19,17,16,15,14,13,12,11,11,10,
        9,9,8,8,7,7,7,6,6,6,5,5,5,4,4,4,3,3,3,3,2,2,2,2,2,1,1,1,1,1,1
    }};
    
    static constexpr PenaltyType default_gap_extension_ = 3;
};

} // namespace octopus

#endif /* hiseq_indel_error_model_hpp */
