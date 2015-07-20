//
//  i_candidate_variant_generator.h
//  Octopus
//
//  Created by Daniel Cooke on 28/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_i_candidate_variant_generator__
#define Octopus_i_candidate_variant_generator__

#include <vector>
#include <cstddef> // std::size_t

#include "common.h"
#include "variant.h"
#include "mappable_set.h"

class AlignedRead;
class GenomicRegion;

class ICandidateVariantGenerator
{
public:
    using RealType = Octopus::ProbabilityType;
    
    // pure virtual functions
    virtual std::vector<Variant> get_candidates(const GenomicRegion&) = 0;
    virtual ~ICandidateVariantGenerator() = default;
    
    virtual void add_read(const AlignedRead&) {};
    virtual void add_reads(std::vector<AlignedRead>::const_iterator first, std::vector<AlignedRead>::const_iterator last) {};
    virtual void add_reads(MappableSet<AlignedRead>::const_iterator first, MappableSet<AlignedRead>::const_iterator last) {};
    virtual void reserve(std::size_t n) {};
    virtual void clear() {};
};

#endif
