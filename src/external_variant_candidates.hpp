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

class ExternalVariantCandidates : public ICandidateVariantGenerator
{
public:
    ExternalVariantCandidates() = delete;
    explicit ExternalVariantCandidates(VcfReader& reader);
    ~ExternalVariantCandidates() override = default;
    
    ExternalVariantCandidates(const ExternalVariantCandidates&)            = default;
    ExternalVariantCandidates& operator=(const ExternalVariantCandidates&) = default;
    ExternalVariantCandidates(ExternalVariantCandidates&&)                 = default;
    ExternalVariantCandidates& operator=(ExternalVariantCandidates&&)      = default;
    
    std::vector<Variant> get_candidates(const GenomicRegion& region) override;
    
private:
    VcfReader& reader_;
};

#endif /* defined(__Octopus__external_variant_candidates__) */
