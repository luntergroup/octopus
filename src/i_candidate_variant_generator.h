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

class AlignedRead;
class GenomicRegion;

class ICandidateVariantGenerator
{
public:
    using RealType     = Octopus::ProbabilityType;
    using ReadIterator = std::vector<AlignedRead>::const_iterator;
    
    // pure virtual functions
    virtual std::vector<Variant> get_candidates(const GenomicRegion&) = 0;
    virtual ~ICandidateVariantGenerator() = default;
    
    virtual void add_read(const AlignedRead&) {};
    virtual void add_reads(ReadIterator first, ReadIterator last) {};
    virtual void reserve(std::size_t n) {};
    virtual void clear() {};
};

#endif
