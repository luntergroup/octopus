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
    VariantFactory(const VariantFactory&)            = delete;
    VariantFactory& operator=(const VariantFactory&) = delete;
    VariantFactory(VariantFactory&&)                 = delete;
    VariantFactory& operator=(VariantFactory&&)      = delete;
    
    static VariantFactory& get_instance()
    {
        static VariantFactory instance;
        return instance;
    }
    
    Variant
    make(std::string contig_name, __uint32_t contig_begin_pos, std::string added, std::string removed) const;
    
private:
    VariantFactory() {};
};

#endif /* defined(__Octopus__variant_factory__) */
