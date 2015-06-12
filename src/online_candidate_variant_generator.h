//
//  online_candidate_variant_generator.h
//  Octopus
//
//  Created by Daniel Cooke on 12/06/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__online_candidate_variant_generator__
#define __Octopus__online_candidate_variant_generator__

#include <vector>

#include "i_candidate_variant_generator.h"

class GenomicRegion;
class Variant;
class AlignedRead;

class OnlineCandidateVariantGenerator : public ICandidateVariantGenerator
{
public:
    OnlineCandidateVariantGenerator() = default;
    ~OnlineCandidateVariantGenerator() override = default;
    
    OnlineCandidateVariantGenerator(const OnlineCandidateVariantGenerator&)            = default;
    OnlineCandidateVariantGenerator& operator=(const OnlineCandidateVariantGenerator&) = default;
    OnlineCandidateVariantGenerator(OnlineCandidateVariantGenerator&&)                 = default;
    OnlineCandidateVariantGenerator& operator=(OnlineCandidateVariantGenerator&&)      = default;
    
    void add_read(const AlignedRead& a_read) override;
    void add_reads(ReadIterator first, ReadIterator last) override;
    std::vector<Variant> get_candidates(const GenomicRegion& a_region) override;
    void reserve(std::size_t n) override;
    void clear() override;
};

#endif /* defined(__Octopus__online_candidate_variant_generator__) */
