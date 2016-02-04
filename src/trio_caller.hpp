//
//  trio_caller.hpp
//  Octopus
//
//  Created by Daniel Cooke on 03/02/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#ifndef trio_caller_hpp
#define trio_caller_hpp

#include <vector>
#include <string>

#include "variant_caller.hpp"
#include "genotype_model.hpp"
#include "trio_genotype_model.hpp"

class GenomicRegion;
class ReadPipe;
class Variant;
class VcfRecord;

namespace Octopus
{
    class TrioVariantCaller : public VariantCaller
    {
    public:
        TrioVariantCaller() = delete;
        
        explicit TrioVariantCaller(const ReferenceGenome& reference,
                                       ReadPipe& read_pipe,
                                       CandidateVariantGenerator&& candidate_generator,
                                       unsigned ploidy,
                                       SampleIdType mother, SampleIdType father,
                                       double min_variant_posterior);
        
        ~TrioVariantCaller() = default;
        
        TrioVariantCaller(const TrioVariantCaller&)            = delete;
        TrioVariantCaller& operator=(const TrioVariantCaller&) = delete;
        TrioVariantCaller(TrioVariantCaller&&)                 = delete;
        TrioVariantCaller& operator=(TrioVariantCaller&&)      = delete;
        
    private:
        const unsigned ploidy_;
        const SampleIdType mother_, father_;
        const double min_variant_posterior_ = 0.95;
        
        std::string do_get_details() const override;
        
        std::vector<VcfRecord> call_variants(const GenomicRegion& region,
                                             const std::vector<Variant>& candidates,
                                             const ReadMap& reads) const override;
    };
    
} // namespace Octopus

#endif /* trio_caller_hpp */
