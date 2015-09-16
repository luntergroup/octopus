//
//  variant_caller.h
//  Octopus
//
//  Created by Daniel Cooke on 26/08/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_variant_caller_h
#define Octopus_variant_caller_h

#include <vector>
#include <iterator>

#include "common.h"
#include "reference_genome.h"
#include "read_manager.h"
#include "read_filter.h"
#include "read_transform.h"
#include "candidate_variant_generator.h"

class GenomicRegion;
class Variant;
class VcfRecord;

class VariantCaller
{
public:
    using ReadFilter = ReadFilter<Octopus::ReadContainer::const_iterator>;
    
    VariantCaller() = delete;
    VariantCaller(ReferenceGenome& reference, ReadManager& read_manager,
                  ReadFilter read_filter, ReadTransform read_transform,
                  CandidateVariantGenerator& candidate_generator);
    virtual ~VariantCaller() = default;
    
    VariantCaller(const VariantCaller&)            = delete;
    VariantCaller& operator=(const VariantCaller&) = delete;
    VariantCaller(VariantCaller&&)                 = delete;
    VariantCaller& operator=(VariantCaller&&)      = delete;
    
    std::vector<VcfRecord> call_variants(const GenomicRegion& region);
    
protected:
    using ReadMap = Octopus::ReadMap;
    
    ReferenceGenome& reference_;
    ReadManager& read_manager_;
    ReadFilter read_filter_;
    ReadTransform read_transform_;
    CandidateVariantGenerator& candidate_generator_;
    
private:
    bool done_calling(const GenomicRegion& region) const noexcept;
    
    virtual GenomicRegion get_init_region(const GenomicRegion& region) = 0;
    virtual GenomicRegion get_next_region(const GenomicRegion& current_region) = 0;
    virtual std::vector<VcfRecord> call_variants(const GenomicRegion& region,
                                                 const std::vector<Variant>& candidates,
                                                 const ReadMap& reads) = 0;
};

#endif
