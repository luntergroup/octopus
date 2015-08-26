//
//  cancer_genotype_model.cpp
//  Octopus
//
//  Created by Daniel Cooke on 26/08/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "cancer_genotype_model.h"

namespace Octopus
{
    // public methods
    
    CancerGenotypeModel::CancerGenotypeModel(unsigned num_samples, SampleIdType normal_sample_id)
    :
    num_samples_ {num_samples},
    normal_sample_id_ {std::move(normal_sample_id)}
    {}
    
    // private methods
    
    CancerGenotypeModel::GenotypeProbabilities
    CancerGenotypeModel::do_evaluate(const std::vector<Haplotype>& haplotypes, const ReadMap& reads)
    {
        GenotypeProbabilities result {};
        
        return result;
    }
    
} // end namespace Octopus
