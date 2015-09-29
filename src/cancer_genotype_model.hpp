//
//  cancer_genotype_model.hpp
//  Octopus
//
//  Created by Daniel Cooke on 26/08/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__cancer_genotype_model__
#define __Octopus__cancer_genotype_model__

#include "genotype_model.hpp"

namespace Octopus
{
    class CancerGenotypeModel : public GenotypeModel
    {
    public:
        CancerGenotypeModel(unsigned num_samples, SampleIdType normal_sample_id,
                            unsigned max_em_iterations = 100, double em_epsilon = 0.001);
        
    private:
        unsigned max_em_iterations_;
        double em_epsilon_;
        
        unsigned num_samples_;
        SampleIdType normal_sample_id_;
        
        GenotypeProbabilities do_evaluate(const std::vector<Haplotype>& haplotypes, const ReadMap& reads) override;
    };
} // namespace Octopus

#endif /* defined(__Octopus__cancer_genotype_model__) */
