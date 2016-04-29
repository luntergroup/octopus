//
//  variant_call.hpp
//  Octopus
//
//  Created by Daniel Cooke on 22/04/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#ifndef variant_call_hpp
#define variant_call_hpp

#include "call.hpp"
#include "allele.hpp"
#include "variant.hpp"

class ReferenceGenome;

namespace Octopus
{
    class VariantCall : public Call
    {
    public:
        using Call::GenotypeCall;
        using Call::PhaseCall;
        
        VariantCall() = delete;
        
        template <typename V, typename T>
        explicit VariantCall(V&& variant, T&& genotype_calls, double quality);
        
        virtual ~VariantCall() = default;
        
        VariantCall(const VariantCall&)            = default;
        VariantCall& operator=(const VariantCall&) = default;
        VariantCall(VariantCall&&)                 = default;
        VariantCall& operator=(VariantCall&&)      = default;
        
        const GenomicRegion& get_region() const noexcept override;
        const Allele& get_reference() const noexcept override;
        
        const Allele& get_alternative() const noexcept;
        
        void replace(const Allele& old, Allele replacement) override;
        void replace_uncalled_genotype_alleles(const Allele& replacement, char ignore) override;
        
        virtual bool parsimonise(char dummy_base) override;
        virtual bool parsimonise(const ReferenceGenome& reference) override;
        
    protected:
        
        Variant variant_;
        
        bool all_genotypes_are_self_contained() const;
        
    private:
        void replace_called_alleles(char old_base, char replacement_base) override;
    };
    
    template <typename V, typename T>
    VariantCall::VariantCall(V&& variant, T&& genotype_calls, double quality)
    : Call {std::forward<T>(genotype_calls), quality}, variant_ {std::forward<V>(variant)}
    {}
} // namespace Octopus

#endif /* variant_call_hpp */
