//
//  variant_benchmarks.cpp
//  Octopus
//
//  Created by Daniel Cooke on 05/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include <iostream>

#include "catch.hpp"
#include "test_common.h"
#include "benchmark_utils.h"

#include "reference_genome.h"
#include "reference_genome_factory.h"
#include "genomic_region.h"
#include "variant.h"
#include "variant_factory.h"
#include "variant_utils.h"

TEST_CASE("variant_benchmark", "[benchmark]")
{
    ReferenceGenomeFactory a_factory {};
    ReferenceGenome human(a_factory.make(human_reference_fasta));
    
    VariantFactory a_variant_factory {};
    
    auto a_region = parse_region("4:3076657-3076660", human);
    
    auto the_sequence = human.get_sequence(a_region);
    
    REQUIRE(the_sequence == "CAG");
    
    auto a_variant = a_variant_factory.make(a_region, the_sequence, std::string {});
    
    auto f_left_align = [&a_variant, &human] () {
        left_align(a_variant, human);
    };
    
    auto time = benchmark<std::chrono::microseconds>(f_left_align, 100).count();
    
    std::cout << "Left alignment time: " << time << "us" << std::endl;
}
