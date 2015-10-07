//
//  cancer_caller.hpp
//  Octopus
//
//  Created by Daniel Cooke on 16/09/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__cancer_caller__
#define __Octopus__cancer_caller__

#include <string>

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
        CancerVariantCaller(ReferenceGenome& reference, CandidateVariantGenerator& candidate_generator,
                            RefCallType refcalls, double min_variant_posterior, double min_somatic_posterior,
                            double min_refcall_posterior, const SampleIdType& normal_sample);
        ~CancerVariantCaller() = default;
        
        CancerVariantCaller(const CancerVariantCaller&)            = delete;
        CancerVariantCaller& operator=(const CancerVariantCaller&) = delete;
        CancerVariantCaller(CancerVariantCaller&&)                 = delete;
        CancerVariantCaller& operator=(CancerVariantCaller&&)      = delete;
        
    private:
        const SampleIdType normal_sample_;
        const double min_variant_posterior_ = 0.95;
        const double min_somatic_mutation_posterior_ = 0.9;
        const double min_refcall_posterior_ = 0.5;
        
        std::string do_get_details() const override;
        
        GenomicRegion get_init_region(const GenomicRegion& region, const ReadMap& reads,
                                      const std::vector<Variant>& candidates) override;
        GenomicRegion get_next_region(const GenomicRegion& current_region, const ReadMap& reads,
                                      const std::vector<Variant>& candidates) override;
        std::vector<VcfRecord> call_variants(const GenomicRegion& region, const std::vector<Variant>& candidates,
                                             const ReadMap& reads) override;
    };
    
} // namespace Octopus

#endif /* defined(__Octopus__cancer_caller__) */
