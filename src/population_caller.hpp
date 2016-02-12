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
#include <string>

#include "variant_caller.hpp"

#include "genotype_model.hpp"
#include "population_genotype_model.hpp"

class GenomicRegion;
class ReadPipe;
class Variant;
class VcfRecord;

namespace Octopus
{
    class PopulationVariantCaller : public VariantCaller
    {
    public:
        PopulationVariantCaller() = delete;
        
        explicit PopulationVariantCaller(const ReferenceGenome& reference,
                                         ReadPipe& read_pipe,
                                         CandidateVariantGenerator&& candidate_generator,
                                         RefCallType refcalls,
                                         double min_variant_posterior,
                                         double min_refcall_posterior,
                                         unsigned ploidy);
        
        ~PopulationVariantCaller() = default;
        
        PopulationVariantCaller(const PopulationVariantCaller&)            = delete;
        PopulationVariantCaller& operator=(const PopulationVariantCaller&) = delete;
        PopulationVariantCaller(PopulationVariantCaller&&)                 = delete;
        PopulationVariantCaller& operator=(PopulationVariantCaller&&)      = delete;
        
    private:
        mutable GenotypeModel::Population genotype_model_;
        
        const unsigned ploidy_;
        const double min_variant_posterior_ = 0.95;
        const double min_refcall_posterior_ = 0.5;
        
        std::vector<VcfRecord> call_variants(const GenomicRegion& region,
                                             const std::vector<Variant>& candidates,
                                             const ReadMap& reads) const override;
    };
    
} // namespace Octopus

#endif /* defined(__Octopus__population_caller__) */
