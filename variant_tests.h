//
//  variant_tests.h
//  Octopus
//
//  Created by Daniel Cooke on 06/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_variant_tests_h
#define Octopus_variant_tests_h

#include <iostream>
#include <assert.h>

#include "genome_region.h"
#include "variant.h"
#include "variant_factory.h"

void allocation_test()
{
    VariantFactory a_variant_factory {};
    
    GenomeRegion ref_region("chr1", 0, 1000);
    
    auto a_variant = a_variant_factory.make(ref_region, "C", "CAT");
    
    std::cout << (*a_variant).get_prior_probability() << std::endl;
}

#endif
