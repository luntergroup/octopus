//
//  variant_factory.h
//  Octopus
//
//  Created by Daniel Cooke on 05/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__variant_factory__
#define __Octopus__variant_factory__

#include "variant.h"

class VariantFactory
{
public:
    VariantFactory() = default;
    
    VariantFactory(const VariantFactory&)            = default;
    VariantFactory& operator=(const VariantFactory&) = default;
    VariantFactory(VariantFactory&&)                 = default;
    VariantFactory& operator=(VariantFactory&&)      = default;
    
    Variant
    make(std::string contig_name, __uint32_t contig_begin_pos, std::string added, std::string removed) const;
    
private:
};

#endif /* defined(__Octopus__variant_factory__) */
