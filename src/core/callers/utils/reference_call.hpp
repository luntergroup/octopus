// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#ifndef reference_call_hpp
#define reference_call_hpp

#include <utility>

#include "call.hpp"

namespace octopus
{
    class ReferenceCall : public Call
    {
    public:
        ReferenceCall() = delete;
        
        template <typename A>
        ReferenceCall(A&& reference, Phred<double> quality);
        
        virtual ~ReferenceCall() = default;
        
        ReferenceCall(const ReferenceCall&)            = default;
        ReferenceCall& operator=(const ReferenceCall&) = default;
        ReferenceCall(ReferenceCall&&)                 = default;
        ReferenceCall& operator=(ReferenceCall&&)      = default;
        
        const GenomicRegion& mapped_region() const noexcept override;
        
        const Allele& reference() const noexcept override;
        
        void replace(const Allele& old, Allele replacement) override;
        void replace_uncalled_genotype_alleles(const Allele& replacement, char ignore) override;
        
        void decorate(VcfRecord::Builder& record) const override;
        
    private:
        Allele reference_;
        
        void replace_called_alleles(char old_base, char replacement_base) override;
    };
    
    template <typename A>
    ReferenceCall::ReferenceCall(A&& reference, Phred<double> quality)
    : Call {quality}, reference_ {std::forward<A>(reference)}
    {}
} // namespace octopus

#endif
