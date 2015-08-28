//
//  haplotype_prior_model.h
//  Octopus
//
//  Created by Daniel Cooke on 26/08/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_haplotype_prior_model_h
#define Octopus_haplotype_prior_model_h

class Haplotype;

class HaplotypePriorModel
{
public:
    double evaluate(const Haplotype& haplotype);
    
    virtual ~HaplotypePriorModel() = default;
    
private:
    virtual double do_evaluate(const Haplotype& haplotype) = 0;
};

inline double HaplotypePriorModel::evaluate(const Haplotype& haplotype)
{
    return do_evaluate(haplotype);
}

#endif
