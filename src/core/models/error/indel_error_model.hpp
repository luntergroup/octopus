// Copyright (c) 2015-2018 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef indel_error_model_hpp
#define indel_error_model_hpp

#include <vector>
#include <array>
#include <cstdint>
#include <memory>

namespace octopus {

class Haplotype;

class IndelErrorModel
{
public:
    using PenaltyType = std::int8_t;
    using PenaltyVector = std::vector<PenaltyType>;
    
    virtual ~IndelErrorModel() = default;
    
    std::unique_ptr<IndelErrorModel> clone() const;
    PenaltyType evaluate(const Haplotype& haplotype, PenaltyVector& gap_open_penalties) const;
    
private:
    virtual std::unique_ptr<IndelErrorModel> do_clone() const = 0;
    virtual PenaltyType do_evaluate(const Haplotype& haplotype, PenaltyVector& gap_open_penalties) const = 0;
};

} // namespace octopus

#endif
