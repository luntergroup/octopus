//
//  reference_genome_implementor.h
//  Octopus
//
//  Created by Daniel Cooke on 10/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_reference_genome_implementor_h
#define Octopus_reference_genome_implementor_h

#include <string>
#include <vector>
#include <cstdint>

#include "genomic_region.h"

class IReferenceGenomeImplementor
{
public:
    virtual std::string get_reference_name() = 0;
    virtual std::vector<std::string> get_contig_names() = 0;
    virtual std::uint_fast32_t get_contig_size(std::string contig_name) = 0;
    virtual std::string get_sequence(const GenomicRegion& a_region) = 0;
    virtual ~IReferenceGenomeImplementor() = 0;
};

inline IReferenceGenomeImplementor::~IReferenceGenomeImplementor() {}

#endif
