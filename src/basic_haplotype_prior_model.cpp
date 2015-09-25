//
//  basic_haplotype_prior_model.cpp
//  Octopus
//
//  Created by Daniel Cooke on 26/08/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "basic_haplotype_prior_model.hpp"

// public methods

BasicBasicHaplotypePriorModel::BasicBasicHaplotypePriorModel(double transition_rate, double transversion_rate)
:
transition_rate_ {transition_rate},
transversion_rate_ {transversion_rate}
{}

// private methods

double BasicBasicHaplotypePriorModel::do_evaluate(const Haplotype& haplotype)
{
    return 1.0;
}
