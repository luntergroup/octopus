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
#include <functional>

#include "i_candidate_variant_generator.hpp"

class GenomicRegion;
class ReferenceGenome;
class Variant;
class AlignedRead;

namespace octopus {
    
class OnlineCandidateVariantGenerator : public ICandidateVariantGenerator
{
public:
    OnlineCandidateVariantGenerator() = delete;
    
    OnlineCandidateVariantGenerator(const ReferenceGenome& reference, Variant::RegionType::Size max_variant_size = 100);
    
    OnlineCandidateVariantGenerator(const OnlineCandidateVariantGenerator&)            = default;
    OnlineCandidateVariantGenerator& operator=(const OnlineCandidateVariantGenerator&) = default;
    OnlineCandidateVariantGenerator(OnlineCandidateVariantGenerator&&)                 = default;
    OnlineCandidateVariantGenerator& operator=(OnlineCandidateVariantGenerator&&)      = default;
    
    ~OnlineCandidateVariantGenerator() override = default;
    
    std::vector<Variant> generate_candidates(const GenomicRegion& region) override;

private:
    std::reference_wrapper<const ReferenceGenome> reference_;
    Variant::RegionType::Size max_variant_size_;
};

} // namespace octopus

#endif /* defined(__Octopus__online_candidate_variant_generator__) */
