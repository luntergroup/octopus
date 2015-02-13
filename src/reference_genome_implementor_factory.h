//
//  reference_genome_implementor_factory.h
//  Octopus
//
//  Created by Daniel Cooke on 10/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__reference_genome_implementor_factory__
#define __Octopus__reference_genome_implementor_factory__

#include <string>
#include <memory>

#include "reference_genome_implementor.h"

class ReferenceGenomeImplementorFactory
{
public:
    ReferenceGenomeImplementorFactory(const ReferenceGenomeImplementorFactory&)            = default;
    ReferenceGenomeImplementorFactory& operator=(const ReferenceGenomeImplementorFactory&) = default;
    ReferenceGenomeImplementorFactory(ReferenceGenomeImplementorFactory&&)                 = default;
    ReferenceGenomeImplementorFactory& operator=(ReferenceGenomeImplementorFactory&&)      = default;
    
    std::unique_ptr<IReferenceGenomeImplementor> make(std::string genome_file_path) const;
};

#endif /* defined(__Octopus__reference_genome_implementor_factory__) */
