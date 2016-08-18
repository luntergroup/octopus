// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef indel_error_model_hpp
#define indel_error_model_hpp

#include <vector>
#include <array>
#include <cstdint>

namespace octopus {

class Haplotype;

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
    static constexpr std::array<PenaltyType, 51> homopolymerErrors_ =
    {{
        60,60,50,45,41,36,31,25,22,20,19,17,16,15,14,13,12,11,11,10,
        9,9,8,8,7,7,7,6,6,6,5,5,5,4,4,4,3,3,3,3,2,2,2,2,2,1,1,1,1,1,1
    }};
    
    static constexpr std::array<PenaltyType, 51> diNucleotideTandemRepeatErrors_ =
    {{
        60,60,48,45,43,41,39,35,31,28,25,21,19,17,15,13,12,11,11,10,
        9,9,8,8,7,7,7,6,6,6,5,5,5,4,4,4,3,3,3,3,2,2,2,2,2,1,1,1,1,1,1
    }};
    
    static constexpr std::array<PenaltyType, 51> triNucleotideTandemRepeatErrors_ =
    {{
        60,60,50,48,46,45,42,39,35,31,28,25,22,20,16,14,13,12,12,11,
        10,9,8,8,7,7,7,6,6,6,5,5,5,4,4,4,3,3,3,3,2,2,2,2,2,1,1,1,1,1,1
    }};
    
    static constexpr std::array<PenaltyType, 51> polyNucleotideTandemRepeatErrors_ =
    {{
        60,60,51,45,45,45,45,45,23,20,19,17,16,15,14,13,12,11,11,10,
        9,9,8,8,7,7,7,6,6,6,5,5,5,4,4,4,3,3,3,3,2,2,2,2,2,1,1,1,1,1,1
    }};
    
    static constexpr PenaltyType defaultGapExtension_ = 3;
};

} // namespace octopus

#endif
