//
//  population_caller.hpp
//  Octopus
//
//  Created by Daniel Cooke on 15/09/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__population_caller__
#define __Octopus__population_caller__

#include <vector>

#include "variant_caller.hpp"

class GenomicRegion;
class ReadManager;
class ReadTransform;
class Variant;
class VcfRecord;

namespace Octopus
{
    class PopulationVariantCaller : public VariantCaller
    {
    public:
        PopulationVariantCaller() = delete;
        PopulationVariantCaller(ReferenceGenome& reference, ReadManager& read_manager,
                                ReadFilter read_filter, ReadTransform read_transform,
                                CandidateVariantGenerator& candidate_generator, unsigned ploidy = 2);
        ~PopulationVariantCaller() = default;
        
        PopulationVariantCaller(const PopulationVariantCaller&)            = delete;
        PopulationVariantCaller& operator=(const PopulationVariantCaller&) = delete;
        PopulationVariantCaller(PopulationVariantCaller&&)                 = delete;
        PopulationVariantCaller& operator=(PopulationVariantCaller&&)      = delete;
        
    private:
        const unsigned ploidy_;
        const double min_posterior_ = 0.95;
        const bool make_ref_calls_ = true;
        
        GenomicRegion get_init_region(const GenomicRegion& region) override;
        GenomicRegion get_next_region(const GenomicRegion& current_region) override;
        std::vector<VcfRecord> call_variants(const GenomicRegion& region, const std::vector<Variant>& candidates,
                                             const ReadMap& reads) override;
    };
    
} // namespace Octopus

#endif /* defined(__Octopus__population_caller__) */
