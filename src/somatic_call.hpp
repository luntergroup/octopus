//
//  somatic_call.hpp
//  Octopus
//
//  Created by Daniel Cooke on 21/04/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#ifndef somatic_call_hpp
#define somatic_call_hpp

#include "variant_call.hpp"

namespace Octopus
{
    class SomaticCall : public VariantCall
    {
    public:
        using VariantCall::GenotypeCall;
        using VariantCall::PhaseCall;
        
        SomaticCall() = delete;
        
        ~SomaticCall() = default;
        
        SomaticCall(const SomaticCall&)            = default;
        SomaticCall& operator=(const SomaticCall&) = default;
        SomaticCall(SomaticCall&&)                 = default;
        SomaticCall& operator=(SomaticCall&&)      = default;
        
        void decorate(VcfRecord::Builder& record) const override;
        
    private:
        
    };
} // namespace Octopus

#endif /* somatic_call_hpp */
