// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

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
    
    GermlineVariantCall(const GermlineVariantCall&)            = default;
    GermlineVariantCall& operator=(const GermlineVariantCall&) = default;
    GermlineVariantCall(GermlineVariantCall&&)                 = default;
    GermlineVariantCall& operator=(GermlineVariantCall&&)      = default;
    
    virtual ~GermlineVariantCall() = default;
    
    void decorate(VcfRecord::Builder& record) const override;
};

template <typename V, typename T>
GermlineVariantCall::GermlineVariantCall(V&& variant, T&& genotype_calls, Phred<double> quality)
: VariantCall {std::forward<V>(variant), std::forward<T>(genotype_calls), quality}
{}

} // namespace octopus

#endif
