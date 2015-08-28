//
//  basic_haplotype_prior_model.h
//  Octopus
//
//  Created by Daniel Cooke on 26/08/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__basic_haplotype_prior_model__
#define __Octopus__basic_haplotype_prior_model__

#include "haplotype_prior_model.h"

class BasicBasicHaplotypePriorModel : public HaplotypePriorModel
{
public:
    BasicBasicHaplotypePriorModel() = default;
    BasicBasicHaplotypePriorModel(double transition_rate, double transversion_rate);
    
private:
    double transition_rate_   = 0.000222;
    double transversion_rate_ = 0.000111;
    
    double do_evaluate(const Haplotype& haplotype) override;
};

#endif /* defined(__Octopus__basic_haplotype_prior_model__) */
