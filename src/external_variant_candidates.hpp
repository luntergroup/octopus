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
#include <memory>

#include <boost/filesystem.hpp>

#include "i_candidate_variant_generator.hpp"
#include "vcf_reader.hpp"

class GenomicRegion;

namespace Octopus {

class ExternalCandidateVariantGenerator : public ICandidateVariantGenerator
{
public:
    ExternalCandidateVariantGenerator() = delete;
    
    ExternalCandidateVariantGenerator(boost::filesystem::path path);
    ExternalCandidateVariantGenerator(std::unique_ptr<const VcfReader> reader);
    ExternalCandidateVariantGenerator(const std::shared_ptr<const VcfReader>& reader);
    
    ~ExternalCandidateVariantGenerator() override = default;
    
    ExternalCandidateVariantGenerator(const ExternalCandidateVariantGenerator&)            = default;
    ExternalCandidateVariantGenerator& operator=(const ExternalCandidateVariantGenerator&) = default;
    ExternalCandidateVariantGenerator(ExternalCandidateVariantGenerator&&)                 = default;
    ExternalCandidateVariantGenerator& operator=(ExternalCandidateVariantGenerator&&)      = default;
    
    std::vector<Variant> generate_candidates(const GenomicRegion& region) override;
    
private:
    std::shared_ptr<const VcfReader> reader_;
};

} // namespace Octopus

#endif /* defined(__Octopus__external_variant_candidates__) */
