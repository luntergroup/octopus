//
//  variant_reader_tests.cpp
//  Octopus
//
//  Created by Daniel Cooke on 15/05/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "catch.hpp"

#include <iostream>
#include <string>

#include "test_common.h"
#include "genomic_region.h"
#include "variant.h"
#include "variant_file_reader.h"
#include "variant_file_factory.h"

using std::cout;
using std::endl;

TEST_CASE("can read", "[variant_file_reader]")
{
    VariantFileFactory a_factory {};
    
    VariantFileReader reader {a_factory.make_reader(sample_vcf)};
    
    auto variants = reader.fetch_variants(GenomicRegion {"X", 10000, 100000});
}
