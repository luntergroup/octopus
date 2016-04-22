//
//  variant_call.cpp
//  Octopus
//
//  Created by Daniel Cooke on 22/04/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#include "variant_call.hpp"

#include "mappable.hpp"

namespace Octopus
{
    const GenomicRegion& VariantCall::get_region() const noexcept
    {
        return mapped_region(variant_);
    }
    
    const Allele& VariantCall::get_reference() const noexcept
    {
        return variant_.get_ref_allele();
    }
    
    const Allele& VariantCall::get_alternative() const noexcept
    {
        return variant_.get_alt_allele();
    }
} // namespace Octopus
