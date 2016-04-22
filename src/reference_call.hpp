//
//  reference_call.hpp
//  Octopus
//
//  Created by Daniel Cooke on 21/04/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#ifndef reference_call_hpp
#define reference_call_hpp

#include <utility>

#include "call.hpp"

namespace Octopus
{
    class ReferenceCall : public Call
    {
    public:
        ReferenceCall() = delete;
        
        template <typename A>
        ReferenceCall(A&& reference, double quality);
        
        ~ReferenceCall() = default;
        
        ReferenceCall(const ReferenceCall&)            = default;
        ReferenceCall& operator=(const ReferenceCall&) = default;
        ReferenceCall(ReferenceCall&&)                 = default;
        ReferenceCall& operator=(ReferenceCall&&)      = default;
        
        const GenomicRegion& get_region() const noexcept override;
        const Allele& get_reference() const noexcept override;
        
        void decorate(VcfRecord::Builder& record) const override;
        
    private:
        Allele reference_;
    };
    
    template <typename A>
    ReferenceCall::ReferenceCall(A&& reference, double quality)
    :
    Call {quality},
    reference_ {std::forward<A>(reference)}
    {}
} // namespace Octopus

#endif /* reference_call_hpp */
