//
//  read_indel_error_model.hpp
//  Octopus
//
//  Created by Daniel Cooke on 17/02/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#ifndef read_indel_error_model_hpp
#define read_indel_error_model_hpp

#include <vector>
#include <array>
#include <cstdint>

#include "haplotype.hpp"

namespace Octopus
{
    class ReadIndelErrorModel
    {
    public:
        using PenaltyType = std::int8_t;
        
        ReadIndelErrorModel() = default;
        ~ReadIndelErrorModel() = default;
        
        ReadIndelErrorModel(const ReadIndelErrorModel&)            = default;
        ReadIndelErrorModel& operator=(const ReadIndelErrorModel&) = default;
        ReadIndelErrorModel(ReadIndelErrorModel&&)                 = default;
        ReadIndelErrorModel& operator=(ReadIndelErrorModel&&)      = default;
        
        PenaltyType calculate_gap_extension_penalty(const Haplotype& haplotype) const noexcept;
        
        void fill_gap_open_penalties(const Haplotype& haplotype, std::vector<PenaltyType>& result) const;
        
    private:
        static constexpr std::array<PenaltyType, 50> Homopolymer_errors_ =
        {{
            45,42,41,39,37,32,28,23,20,19,17,16,15,14,13,12,11,11,10,
            9,9,8,8,7,7,7,6,6,6,5,5,5,4,4,4,3,3,3,3,2,2,2,2,2,1,1,1,1,1,1
        }};
        
        static constexpr PenaltyType gap_extension_ = 3;
    };
} // namespace Octopus

#endif /* read_indel_error_model_hpp */
