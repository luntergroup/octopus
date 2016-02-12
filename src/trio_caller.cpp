//
//  trio_caller.cpp
//  Octopus
//
//  Created by Daniel Cooke on 03/02/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#include "trio_caller.hpp"

#include <utility>

#include "read_pipe.hpp"
#include "vcf_record.hpp"

namespace Octopus
{
    // public methods
    
    TrioVariantCaller::TrioVariantCaller(const ReferenceGenome& reference,
                                         ReadPipe& read_pipe,
                                         CandidateVariantGenerator&& candidate_generator,
                                         unsigned ploidy,
                                         SampleIdType mother, SampleIdType father,
                                         double min_variant_posterior)
    :
    VariantCaller {reference, read_pipe, std::move(candidate_generator), RefCallType::None},
    ploidy_ {ploidy},
    mother_ {std::move(mother)},
    father_ {std::move(father)},
    min_variant_posterior_ {min_variant_posterior}
    {}
    
    // private methods
    
    namespace {
        using GM = GenotypeModel::Trio;
        
    } // namespace
    
    std::vector<VcfRecord>
    TrioVariantCaller::call_variants(const GenomicRegion& region,
                                     const std::vector<Variant>& candidates,
                                     const ReadMap& reads) const
    {
        std::vector<VcfRecord> result {};
        return result;
    }
    
} // namespace Octopus
