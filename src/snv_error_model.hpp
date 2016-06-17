//
//  snv_error_model.hpp
//  Octopus
//
//  Created by Daniel Cooke on 15/06/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#ifndef snv_error_model_hpp
#define snv_error_model_hpp

#include <vector>
#include <array>
#include <cstdint>

class Haplotype;

namespace Octopus
{
    class SnvErrorModel
    {
    public:
        using PenaltyType = std::int8_t;
        
        using PenaltyVector = std::vector<PenaltyType>;
        
        SnvErrorModel() = default;
        
        virtual ~SnvErrorModel() = default;
        
        void evaluate(const Haplotype& haplotype,
                      PenaltyVector& forward_snv_priors,
                      PenaltyVector& reverse_snv_priors) const;
        
    private:
        static constexpr std::array<PenaltyType, 51> Homopolymer_errors_ =
        {{
            125,125,60,55,40,30,25,20,16,14,13,12,12,11,11,10,10,10,9,8,
            7,7,7,6,6,6,6,6,6,5,5,5,5,4,4,4,3,3,3,3,2,2,2,2,2,1,1,1,1,1,1
        }};
    };
} // namespace Octopus

#endif /* snv_error_model_hpp */
