//
//  variant_factory.h
//  Octopus
//
//  Created by Daniel Cooke on 05/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__variant_factory__
#define __Octopus__variant_factory__

#include <memory>

#include "variant.h"

class VariantFactory
{
public:
    VariantFactory() = delete;
    
    std::unique_ptr<Variant>
    make(std::string sequence_name, __uint32_t sequence_start_pos, std::string added, std::string removed) const;
    
private:
};

#endif /* defined(__Octopus__variant_factory__) */
