//
//  pedigree_caller.cpp
//  Octopus
//
//  Created by Daniel Cooke on 30/10/2015.
//  Copyright © 2015 Oxford University. All rights reserved.
//

#include "pedigree_caller.hpp"

#include <utility>

#include "read_pipe.hpp"
#include "vcf_record.hpp"

namespace octopus
{
    // public methods
    
//    PedigreeVariantCaller::PedigreeVariantCaller(const ReferenceGenome& reference,
//                                                 ReadPipe& read_pipe,
//                                                 CandidateVariantGenerator&& candidate_generator,
//                                                 unsigned ploidy,
//                                                 SampleName mother, SampleName father,
//                                                 double min_variant_posterior)
//    :
//    //VariantCaller {reference, read_pipe, std::move(candidate_generator), RefCallType::None},
//    ploidy_ {ploidy},
//    mother_ {std::move(mother)},
//    father_ {std::move(father)},
//    min_variant_posterior_ {min_variant_posterior}
//    {}
    
    // private methods
    
    namespace
    {
        using GM = model::Pedigree;
    } // namespace

} // namespace octopus
