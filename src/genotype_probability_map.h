//
//  genotype_probability_map.h
//  Octopus
//
//  Created by Daniel Cooke on 04/02/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#ifndef genotype_probability_map_h
#define genotype_probability_map_h

#include "common.hpp"
#include "matrix_map.hpp"

template <typename GenotypeTp>
using GenotypeProbabilityMap = MatrixMap<SampleIdType, GenotypeTp, double>;

#endif /* genotype_probability_map_h */
