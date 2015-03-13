//
//  variant_candidate_generator_interface.h
//  Octopus
//
//  Created by Daniel Cooke on 28/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_i_variant_candidate_generator_h
#define Octopus_i_variant_candidate_generator_h

#include <vector>
#include <cstddef> // std::size_t

class Variant;
class AlignedRead;
class GenomicRegion;

class IVariantCandidateGenerator
{
public:
    using ReadIterator = std::vector<AlignedRead>::const_iterator;
    
    virtual std::vector<Variant> get_candidates(const GenomicRegion&) = 0;
    virtual void add_read(const AlignedRead&) = 0;
    virtual void add_reads(ReadIterator first, ReadIterator last) = 0;
    virtual void reserve(std::size_t n) = 0;
    virtual void clear() = 0;
    virtual double get_variant_detection_probability(const Variant& a_variant) = 0;
    virtual ~IVariantCandidateGenerator() = default;
};

#endif
