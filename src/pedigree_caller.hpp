//
//  Pedigree_caller.hpp
//  Octopus
//
//  Created by Daniel Cooke on 30/10/2015.
//  Copyright © 2015 Oxford University. All rights reserved.
//

#ifndef pedigree_caller_hpp
#define pedigree_caller_hpp

#include <vector>
#include <string>

#include "variant_caller.hpp"

class GenomicRegion;
class ReadManager;
class ReadTransform;
class Variant;
class VcfRecord;

namespace Octopus
{
    class PedigreeVariantCaller : public VariantCaller
    {
    public:
        PedigreeVariantCaller() = delete;
        
        explicit PedigreeVariantCaller(const ReferenceGenome& reference,
                                       CandidateVariantGenerator&& candidate_generator,
                                       unsigned ploidy,
                                       SampleIdType mother, SampleIdType father,
                                       double min_variant_posterior);
        
        ~PedigreeVariantCaller() = default;
        
        PedigreeVariantCaller(const PedigreeVariantCaller&)            = delete;
        PedigreeVariantCaller& operator=(const PedigreeVariantCaller&) = delete;
        PedigreeVariantCaller(PedigreeVariantCaller&&)                 = delete;
        PedigreeVariantCaller& operator=(PedigreeVariantCaller&&)      = delete;
        
    private:
        const unsigned ploidy_;
        const SampleIdType mother_, father_;
        const double min_variant_posterior_ = 0.95;
        
        std::string do_get_details() const override;
        
        std::vector<VcfRecord> call_variants(const GenomicRegion& region, const std::vector<Variant>& candidates,
                                             const ReadMap& reads) override;
    };
    
} // namespace Octopus

#endif /* pedigree_caller_hpp */
