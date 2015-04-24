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

#include "variant.h"

class AlignedRead;
class GenomicRegion;

class ICandidateVariantGenerator
{
public:
    using RealType     = Variant::RealType;
    using ReadIterator = std::vector<AlignedRead>::const_iterator;
    
    virtual std::vector<Variant> get_candidates(const GenomicRegion&) = 0;
    virtual void add_read(const AlignedRead&) = 0;
    virtual void add_reads(ReadIterator first, ReadIterator last) = 0;
    virtual void reserve(std::size_t n) = 0;
    virtual void clear() = 0;
    virtual ~ICandidateVariantGenerator() = default;
};

#endif
