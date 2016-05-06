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

#include "allele.hpp"

namespace Octopus
{
    class GermlineVariantCall : public VariantCall
    {
    public:
        using VariantCall::GenotypeCall;
        using VariantCall::PhaseCall;
        
        GermlineVariantCall() = delete;
        
        template <typename V, typename T>
        explicit GermlineVariantCall(V&& variant, T&& genotype_calls, double quality);
        
        virtual ~GermlineVariantCall() = default;
        
        GermlineVariantCall(const GermlineVariantCall&)            = default;
        GermlineVariantCall& operator=(const GermlineVariantCall&) = default;
        GermlineVariantCall(GermlineVariantCall&&)                 = default;
        GermlineVariantCall& operator=(GermlineVariantCall&&)      = default;
        
        void decorate(VcfRecord::Builder& record) const override;
    };
    
    template <typename V, typename T>
    GermlineVariantCall::GermlineVariantCall(V&& variant, T&& genotype_calls, double quality)
    : VariantCall {std::forward<V>(variant), std::forward<T>(genotype_calls), quality}
    {}
} // namespace Octopus

#endif /* germline_variant_call_hpp */
