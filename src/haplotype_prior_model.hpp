//
//  haplotype_prior_model.hpp
//  Octopus
//
//  Created by Daniel Cooke on 26/08/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_haplotype_prior_model_hpp
#define Octopus_haplotype_prior_model_hpp

class Haplotype;

class HaplotypePriorModel
{
public:
    HaplotypePriorModel() = delete;
    
    double evaluate(const Haplotype& haplotype);
};

inline double HaplotypePriorModel::evaluate(const Haplotype& haplotype)
{
    return 0;
}

#endif
