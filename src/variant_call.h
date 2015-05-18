//
//  variant_call.h
//  Octopus
//
//  Created by Daniel Cooke on 18/05/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_variant_call_h
#define Octopus_variant_call_h

#include "common.h"
#include "genomic_region.h"
#include "comparable.h"
#include "mappable.h"
#include "variant.h"

class VariantCall : public Comparable<VariantCall>, public Mappable<VariantCall>
{
public:
    using RealType = Octopus::ProbabilityType;
    
    VariantCall() = delete;
    ~VariantCall() = default;
    
    VariantCall(const VariantCall&)            = default;
    VariantCall& operator=(const VariantCall&) = default;
    VariantCall(VariantCall&&)                 = default;
    VariantCall& operator=(VariantCall&&)      = default;
    
    
private:
    Variant the_variant_;
    RealType the_call_correctness_probability_;
}

#endif
