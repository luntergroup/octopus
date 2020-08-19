// Copyright (c) 2015-2020 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef snv_error_model_hpp
#define snv_error_model_hpp

#include <vector>
#include <cstdint>
#include <memory>

namespace octopus {

class Haplotype;

class SnvErrorModel
{
public:
    using MutationVector = std::vector<char>;
    using PenaltyType    = std::int8_t;
    using PenaltyVector  = std::vector<PenaltyType>;
    
    virtual ~SnvErrorModel() = default;
    
    std::unique_ptr<SnvErrorModel> clone() const;
    void evaluate(const Haplotype& haplotype,
                  MutationVector& forward_snv_mask, PenaltyVector& forward_snv_priors,
                  MutationVector& reverse_snv_mask, PenaltyVector& reverse_snv_priors) const;

private:
    virtual std::unique_ptr<SnvErrorModel> do_clone() const = 0;
    virtual void do_evaluate(const Haplotype& haplotype,
                             MutationVector& forward_snv_mask, PenaltyVector& forward_snv_priors,
                             MutationVector& reverse_snv_mask, PenaltyVector& reverse_snv_priors) const = 0;
};

} // namespace octopus

#endif
