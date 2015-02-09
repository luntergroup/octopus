//
//  variant_tests.cpp
//  Octopus
//
//  Created by Daniel Cooke on 09/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include <iostream>

#define CATCH_CONFIG_MAIN
#include "catch.hpp"

#include "genomic_region.h"
#include "variant.h"
#include "variant_factory.h"

TEST_CASE("creation_test", "snp")
{
    VariantFactory a_variant_factory {};
    
    GenomicRegion ref_region("chr1", 0, 1000);
    
    auto a_variant = a_variant_factory.make(ref_region, "C", "CAT");
    
    std::cout << a_variant->get_prior_probability() << std::endl;
    
    REQUIRE(a_variant->get_prior_probability() == 1e-7);
}
