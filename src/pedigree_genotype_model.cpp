//
//  pedigree_genotype_model.cpp
//  Octopus
//
//  Created by Daniel Cooke on 23/11/2015.
//  Copyright © 2015 Oxford University. All rights reserved.
//

#include "pedigree_genotype_model.hpp"

namespace Octopus
{
    namespace GenotypeModel
    {
        namespace {
            
        } // namespace
        
        Pedigree::Latents
        Pedigree::infer_latents(const std::vector<Haplotype>& haplotypes, const ReadMap& reads,
                                const ReferenceGenome& reference)
        {
            Pedigree::Latents result {};
            return result;
        }
    } // namespace GenotypeModel
} // namespace Octopus
