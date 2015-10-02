//
//  cancer_caller.hpp
//  Octopus
//
//  Created by Daniel Cooke on 16/09/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__cancer_caller__
#define __Octopus__cancer_caller__

#include "variant_caller.hpp"

class GenomicRegion;
class ReadManager;
class ReadTransform;
class Variant;
class VcfRecord;

namespace Octopus
{

class CancerVariantCaller : public VariantCaller
{
public:
    CancerVariantCaller() = delete;
    CancerVariantCaller(ReferenceGenome& reference, ReadManager& read_manager,
                        ReadFilter read_filter, ReadTransform read_transform,
                        CandidateVariantGenerator& candidate_generator);
    ~CancerVariantCaller() = default;
    
    CancerVariantCaller(const CancerVariantCaller&)            = delete;
    CancerVariantCaller& operator=(const CancerVariantCaller&) = delete;
    CancerVariantCaller(CancerVariantCaller&&)                 = delete;
    CancerVariantCaller& operator=(CancerVariantCaller&&)      = delete;
    
private:
    const double min_posterior_ = 0.95;
    
    GenomicRegion get_init_region(const GenomicRegion& region) override;
    GenomicRegion get_next_region(const GenomicRegion& current_region) override;
    std::vector<VcfRecord> call_variants(const GenomicRegion& region, const std::vector<Variant>& candidates,
                                         const ReadMap& reads) override;
};
    
} // namespace Octopus

#endif /* defined(__Octopus__cancer_caller__) */
