//
//  indel_error_model.hpp
//  Octopus
//
//  Created by Daniel Cooke on 17/02/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#ifndef indel_error_model_hpp
#define indel_error_model_hpp

#include <vector>
#include <array>
#include <cstdint>

class Haplotype;

namespace Octopus
{
    class IndelErrorModel
    {
    public:
        using PenaltyType = std::int8_t;
        
        using PenaltyVector = std::vector<PenaltyType>;
        
        IndelErrorModel() = default;
        
        IndelErrorModel(const IndelErrorModel&)            = default;
        IndelErrorModel& operator=(const IndelErrorModel&) = default;
        IndelErrorModel(IndelErrorModel&&)                 = default;
        IndelErrorModel& operator=(IndelErrorModel&&)      = default;
        
        virtual ~IndelErrorModel() = default;
        
        PenaltyType evaluate(const Haplotype& haplotype, PenaltyVector& gap_open_penalties) const;
        
    private:
        static constexpr std::array<PenaltyType, 51> Homopolymer_errors_ =
        {{
            60,60,50,45,41,36,31,25,22,20,19,17,16,15,14,13,12,11,11,10,
            9,9,8,8,7,7,7,6,6,6,5,5,5,4,4,4,3,3,3,3,2,2,2,2,2,1,1,1,1,1,1
        }};
        
        static constexpr std::array<PenaltyType, 51> Di_nucleotide_tandem_repeat_errors_ =
        {{
            60,60,48,45,43,41,39,35,31,28,25,21,19,17,15,13,12,11,11,10,
            9,9,8,8,7,7,7,6,6,6,5,5,5,4,4,4,3,3,3,3,2,2,2,2,2,1,1,1,1,1,1
        }};
        
        static constexpr std::array<PenaltyType, 51> Tri_nucleotide_tandem_repeat_errors_ =
        {{
            60,60,50,48,46,45,42,39,35,31,28,25,22,20,16,14,13,12,12,11,
            10,9,8,8,7,7,7,6,6,6,5,5,5,4,4,4,3,3,3,3,2,2,2,2,2,1,1,1,1,1,1
        }};
        
        static constexpr std::array<PenaltyType, 51> Poly_nucleotide_tandem_repeat_errors_ =
        {{
            60,60,51,45,45,45,45,45,23,20,19,17,16,15,14,13,12,11,11,10,
            9,9,8,8,7,7,7,6,6,6,5,5,5,4,4,4,3,3,3,3,2,2,2,2,2,1,1,1,1,1,1
        }};
        
        static constexpr PenaltyType default_gap_extension_ = 3;
    };
} // namespace Octopus

#endif /* indel_error_model_hpp */
