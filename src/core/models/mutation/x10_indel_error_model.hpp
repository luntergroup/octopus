// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef x10_indel_error_model_hpp
#define x10_indel_error_model_hpp

#include "indel_error_model.hpp"

namespace octopus {

class Haplotype;

class X10IndelErrorModel : public IndelErrorModel
{
public:
    using IndelErrorModel::PenaltyType;
    using IndelErrorModel::PenaltyVector;
    
    X10IndelErrorModel() = default;
    
    X10IndelErrorModel(const X10IndelErrorModel&)            = default;
    X10IndelErrorModel& operator=(const X10IndelErrorModel&) = default;
    X10IndelErrorModel(X10IndelErrorModel&&)                 = default;
    X10IndelErrorModel& operator=(X10IndelErrorModel&&)      = default;

private:
    static constexpr std::array<PenaltyType, 50> homopolymerErrors_ =
    {{
     45,42,41,39,37,32,28,23,20,18,16,15,14,13,12,11,10,10,9,
     9,8,8,8,7,7,7,6,6,6,5,5,5,4,4,4,3,3,3,3,2,2,2,2,2,1,1,1,1,1,1
     }};
    static constexpr std::array<PenaltyType, 50> diNucleotideTandemRepeatErrors_ =
    {{
     45,45,45,45,36,31,27,22,19,18,16,15,14,13,12,11,10,10,10,
     9,9,8,8,7,7,7,6,6,6,5,5,5,4,4,4,3,3,3,3,2,2,2,2,2,1,1,1,1,1,1
     }};
    static constexpr std::array<PenaltyType, 50> triNucleotideTandemRepeatErrors_ =
    {{
     45,45,45,45,45,45,28,23,20,19,17,16,15,14,13,12,11,11,10,
     9,9,8,8,7,7,7,6,6,6,5,5,5,4,4,4,3,3,3,3,2,2,2,2,2,1,1,1,1,1,1
     }};
    static constexpr std::array<PenaltyType, 50> polyNucleotideTandemRepeatErrors_ =
    {{
     45,45,45,45,45,45,45,23,20,19,17,16,15,14,13,12,11,11,10,
     9,9,8,8,7,7,7,6,6,6,5,5,5,4,4,4,3,3,3,3,2,2,2,2,2,1,1,1,1,1,1
     }};
    
    static constexpr PenaltyType defaultGapExtension_ = 2;
    
    virtual PenaltyType do_evaluate(const Haplotype& haplotype, PenaltyVector& gap_open_penalties) const;
};
    
} // namespace octopus

#endif
