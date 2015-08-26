//
//  cancer_genotype_model.h
//  Octopus
//
//  Created by Daniel Cooke on 26/08/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__cancer_genotype_model__
#define __Octopus__cancer_genotype_model__

#include "genotype_model.h"

namespace Octopus
{
    class CancerGenotypeModel : public GenotypeModel
    {
    public:
        CancerGenotypeModel(unsigned num_samples, SampleIdType normal_sample_id);
        
    private:
        unsigned num_samples_;
        SampleIdType normal_sample_id_;
        
        GenotypeProbabilities do_evaluate(const std::vector<Haplotype>& haplotypes, const ReadMap& reads) override;
    };
} // end namespace Octopus

#endif /* defined(__Octopus__cancer_genotype_model__) */
