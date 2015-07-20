//
//  external_variant_candidates.h
//  Octopus
//
//  Created by Daniel Cooke on 28/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__external_variant_candidates__
#define __Octopus__external_variant_candidates__

#include <cstddef> // std::size_t

#include "i_candidate_variant_generator.h"
#include "variant_file_reader.h"

class GenomicRegion;

class ExternalVariantCandidates : public ICandidateVariantGenerator
{
public:
    ExternalVariantCandidates() = delete;
    explicit ExternalVariantCandidates(VariantFileReader& a_variant_source, double generator_confidence);
    ~ExternalVariantCandidates() override = default;
    
    ExternalVariantCandidates(const ExternalVariantCandidates&)            = default;
    ExternalVariantCandidates& operator=(const ExternalVariantCandidates&) = default;
    ExternalVariantCandidates(ExternalVariantCandidates&&)                 = default;
    ExternalVariantCandidates& operator=(ExternalVariantCandidates&&)      = default;
    
    std::vector<Variant> get_candidates(const GenomicRegion& a_region) override;
    
private:
    VariantFileReader& a_variant_file_;
    double generator_confidence_;
};

#endif /* defined(__Octopus__external_variant_candidates__) */
