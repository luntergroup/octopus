//
//  reference_genome_implementor_factory.h
//  Octopus
//
//  Created by Daniel Cooke on 10/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__reference_genome_impl_factory__
#define __Octopus__reference_genome_impl_factory__

#include <string>
#include <memory> // std::unique_ptr, std::make_unique

class IReferenceGenomeImpl;

class ReferenceGenomeImplFactory
{
public:
    ReferenceGenomeImplFactory(const ReferenceGenomeImplFactory&)            = default;
    ReferenceGenomeImplFactory& operator=(const ReferenceGenomeImplFactory&) = default;
    ReferenceGenomeImplFactory(ReferenceGenomeImplFactory&&)                 = default;
    ReferenceGenomeImplFactory& operator=(ReferenceGenomeImplFactory&&)      = default;
    
    std::unique_ptr<IReferenceGenomeImpl> make(std::string genome_file_path) const;
};

#endif /* defined(__Octopus__reference_genome_implementor_factory__) */
