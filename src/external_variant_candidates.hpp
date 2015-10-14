//
//  external_variant_candidates.hpp
//  Octopus
//
//  Created by Daniel Cooke on 28/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__external_variant_candidates__
#define __Octopus__external_variant_candidates__

#include <vector>

#include "i_candidate_variant_generator.hpp"
#include "vcf_reader.hpp"

class GenomicRegion;

namespace Octopus {

class ExternalCandidateVariantGenerator : public ICandidateVariantGenerator
{
public:
    ExternalCandidateVariantGenerator() = delete;
    explicit ExternalCandidateVariantGenerator(VcfReader&& reader);
    ~ExternalCandidateVariantGenerator() override = default;
    
    ExternalCandidateVariantGenerator(const ExternalCandidateVariantGenerator&)            = default;
    ExternalCandidateVariantGenerator& operator=(const ExternalCandidateVariantGenerator&) = default;
    ExternalCandidateVariantGenerator(ExternalCandidateVariantGenerator&&)                 = default;
    ExternalCandidateVariantGenerator& operator=(ExternalCandidateVariantGenerator&&)      = default;
    
    std::vector<Variant> get_candidates(const GenomicRegion& region) override;
    
private:
    VcfReader reader_;
};

} // namespace Octopus

#endif /* defined(__Octopus__external_variant_candidates__) */
