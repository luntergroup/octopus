//
//  germline_variant_call.hpp
//  Octopus
//
//  Created by Daniel Cooke on 21/04/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#ifndef germline_variant_call_hpp
#define germline_variant_call_hpp

#include "variant_call.hpp"

#include <vector>
#include <utility>

#include <core/types/allele.hpp>

namespace octopus {

class GermlineVariantCall : public VariantCall
{
public:
    using VariantCall::GenotypeCall;
    using VariantCall::PhaseCall;
    
    GermlineVariantCall() = delete;
    
    template <typename V, typename T>
    GermlineVariantCall(V&& variant, T&& genotype_calls, Phred<double> quality);
    
    virtual ~GermlineVariantCall() = default;
    
    GermlineVariantCall(const GermlineVariantCall&)            = default;
    GermlineVariantCall& operator=(const GermlineVariantCall&) = default;
    GermlineVariantCall(GermlineVariantCall&&)                 = default;
    GermlineVariantCall& operator=(GermlineVariantCall&&)      = default;
    
    void decorate(VcfRecord::Builder& record) const override;
};

template <typename V, typename T>
GermlineVariantCall::GermlineVariantCall(V&& variant, T&& genotype_calls, Phred<double> quality)
: VariantCall {std::forward<V>(variant), std::forward<T>(genotype_calls), quality}
{}

} // namespace octopus

#endif /* germline_variant_call_hpp */
