//
//  basic_caller.h
//  Octopus
//
//  Created by Daniel Cooke on 15/09/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__basic_caller__
#define __Octopus__basic_caller__

#include <vector>

#include "variant_caller.h"

class GenomicRegion;
class ReadManager;
class ReadTransform;
class Variant;
class VcfRecord;

class BasicVariantCaller : public VariantCaller
{
public:
    BasicVariantCaller() = delete;
    BasicVariantCaller(ReferenceGenome& reference, ReadManager& read_manager,
                       ReadFilter read_filter, ReadTransform read_transform,
                       CandidateVariantGenerator& candidate_generator);
    ~BasicVariantCaller() = default;
    
    BasicVariantCaller(const BasicVariantCaller&)            = delete;
    BasicVariantCaller& operator=(const BasicVariantCaller&) = delete;
    BasicVariantCaller(BasicVariantCaller&&)                 = delete;
    BasicVariantCaller& operator=(BasicVariantCaller&&)      = delete;
    
private:
    GenomicRegion get_init_region(const GenomicRegion& region) override;
    GenomicRegion get_next_region(const GenomicRegion& current_region) override;
    std::vector<VcfRecord> call_variants(const GenomicRegion& region, const std::vector<Variant>& candidates,
                                         const ReadMap& reads) override;
};

#endif /* defined(__Octopus__basic_caller__) */
