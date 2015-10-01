//
//  online_candidate_variant_generator.hpp
//  Octopus
//
//  Created by Daniel Cooke on 12/06/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__online_candidate_variant_generator__
#define __Octopus__online_candidate_variant_generator__

#include <vector>

#include "i_candidate_variant_generator.hpp"

class GenomicRegion;
class ReferenceGenome;
class Variant;
class AlignedRead;

class OnlineCandidateVariantGenerator : public ICandidateVariantGenerator
{
public:
    using SizeType = GenomicRegion::SizeType;
    
    OnlineCandidateVariantGenerator() = delete;
    explicit OnlineCandidateVariantGenerator(ReferenceGenome& reference, SizeType max_variant_size = 100);
    ~OnlineCandidateVariantGenerator() override = default;
    
    OnlineCandidateVariantGenerator(const OnlineCandidateVariantGenerator&)            = default;
    OnlineCandidateVariantGenerator& operator=(const OnlineCandidateVariantGenerator&) = default;
    OnlineCandidateVariantGenerator(OnlineCandidateVariantGenerator&&)                 = default;
    OnlineCandidateVariantGenerator& operator=(OnlineCandidateVariantGenerator&&)      = default;
    
    std::vector<Variant> get_candidates(const GenomicRegion& region) override;

private:
    ReferenceGenome& reference_;
    SizeType max_variant_size_;
};

#endif /* defined(__Octopus__online_candidate_variant_generator__) */
