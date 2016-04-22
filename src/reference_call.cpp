//
//  reference_call.cpp
//  Octopus
//
//  Created by Daniel Cooke on 21/04/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#include "reference_call.hpp"

#include "mappable.hpp"

namespace Octopus
{
    const GenomicRegion& ReferenceCall::get_region() const noexcept
    {
        return mapped_region(reference_);
    }
    
    const Allele& ReferenceCall::get_reference() const noexcept
    {
        return reference_;
    }
    
    void ReferenceCall::decorate(VcfRecord::Builder& record) const
    {
        
    }
} // namespace Octopus
