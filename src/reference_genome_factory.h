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

#include "reference_genome_impl.h"

class ReferenceGenomeFactory
{
public:
    ReferenceGenomeFactory(const ReferenceGenomeFactory&)            = default;
    ReferenceGenomeFactory& operator=(const ReferenceGenomeFactory&) = default;
    ReferenceGenomeFactory(ReferenceGenomeFactory&&)                 = default;
    ReferenceGenomeFactory& operator=(ReferenceGenomeFactory&&)      = default;
    
    std::unique_ptr<IReferenceGenomeImpl> make(std::string file_path,
                                               IReferenceGenomeImpl::SizeType max_cache_size = 0) const;
    
    std::unique_ptr<IReferenceGenomeImpl> make(std::string file_path, std::string index_path,
                                               IReferenceGenomeImpl::SizeType max_cache_size = 0) const;
};

#endif /* defined(__Octopus__reference_genome_implementor_factory__) */
