//
//  trio_genotype_model.cpp
//  Octopus
//
//  Created by Daniel Cooke on 29/01/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#include "trio_genotype_model.hpp"

namespace Octopus
{
    namespace GenotypeModel
    {
        namespace {
            
        } // namespace
        
        Trio::Latents
        Trio::infer_latents(const std::vector<Haplotype>& haplotypes, const ReadMap& reads,
                            const ReferenceGenome& reference)
        {
            Trio::Latents result {};
            return result;
        }
    } // namespace GenotypeModel
} // namespace Octopus
