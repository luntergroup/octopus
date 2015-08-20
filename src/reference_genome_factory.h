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
#include <cstddef>
#include <boost/filesystem/path.hpp>

#include "i_reference_genome_impl.h"
#include "reference_genome.h"

namespace fs = boost::filesystem;

ReferenceGenome make_reference(fs::path file_path, std::size_t max_base_pair_cache = 0, bool is_threaded = false);
//ReferenceGenome make_reference(fs::path file_path, bool is_threaded = false);

class ReferenceGenomeFactory
{
public:
    ReferenceGenomeFactory(const ReferenceGenomeFactory&)            = default;
    ReferenceGenomeFactory& operator=(const ReferenceGenomeFactory&) = default;
    ReferenceGenomeFactory(ReferenceGenomeFactory&&)                 = default;
    ReferenceGenomeFactory& operator=(ReferenceGenomeFactory&&)      = default;
    
    std::unique_ptr<IReferenceGenomeImpl> make(fs::path file_path,
                                               IReferenceGenomeImpl::SizeType max_cache_size = 0) const;
    
    std::unique_ptr<IReferenceGenomeImpl> make(fs::path file_path, std::string index_path,
                                               IReferenceGenomeImpl::SizeType max_cache_size = 0) const;
};

#endif /* defined(__Octopus__reference_genome_implementor_factory__) */
