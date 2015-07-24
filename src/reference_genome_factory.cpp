//
//  reference_genome_implementor_factory.cpp
//  Octopus
//
//  Created by Daniel Cooke on 10/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "reference_genome_factory.h"

#include "fasta.h"
#include "caching_fasta.h"

std::unique_ptr<IReferenceGenomeImpl>
ReferenceGenomeFactory::make(std::string file_path, IReferenceGenomeImpl::SizeType max_cache_size) const
{
    if (max_cache_size == 0) {
        return std::make_unique<Fasta>(std::move(file_path));
    } else {
        return std::make_unique<CachingFasta>(std::move(file_path), max_cache_size);
    }
}

std::unique_ptr<IReferenceGenomeImpl>
ReferenceGenomeFactory::make(std::string file_path, std::string index_path,
                             IReferenceGenomeImpl::SizeType max_cache_size) const
{
    if (max_cache_size == 0) {
        return std::make_unique<Fasta>(std::move(file_path), std::move(index_path));
    } else {
        return std::make_unique<CachingFasta>(std::move(file_path), std::move(index_path), max_cache_size);
    }
}
