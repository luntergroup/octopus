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

//TEST_CASE("variant_benchmark", "[benchmark]")
//{
//    ReferenceGenomeFactory a_factory {};
//    ReferenceGenome human(a_factory.make(human_reference_fasta));
//    
//    VariantFactory a_variant_factory {};
//    
//    auto a_region = parse_region("4:3076657-3076660", human);
//    
//    auto the_sequence = human.get_sequence(a_region);
//    
//    REQUIRE(the_sequence == "CAG");
//    
//    auto a_variant = a_variant_factory.make(a_region, the_sequence, std::string {});
//    
//    auto f_left_align_1 = [&a_variant, &human] () {
//        left_align(a_variant, human, 1);
//    };
//    
//    auto f_left_align_10 = [&a_variant, &human] () {
//        left_align(a_variant, human, 10);
//    };
//    
//    auto f_left_align_20 = [&a_variant, &human] () {
//        left_align(a_variant, human, 20);
//    };
//    
//    auto f_left_align_30 = [&a_variant, &human] () {
//        left_align(a_variant, human, 30);
//    };
//    
//    auto f_left_align_100 = [&a_variant, &human] () {
//        left_align(a_variant, human, 100);
//    };
//    
//    auto extension_1_time = benchmark<std::chrono::microseconds>(f_left_align_1, 100).count();
//    auto extension_10_time = benchmark<std::chrono::microseconds>(f_left_align_10, 100).count();
//    auto extension_20_time = benchmark<std::chrono::microseconds>(f_left_align_20, 100).count();
//    auto extension_30_time = benchmark<std::chrono::microseconds>(f_left_align_30, 100).count();
//    auto extension_100_time = benchmark<std::chrono::microseconds>(f_left_align_100, 100).count();
//    
//    
//    std::cout << "extension_1_time left alignment time: " << extension_1_time << "us" << std::endl;
//    std::cout << "extension_10_time left alignment time: " << extension_10_time << "us" << std::endl;
//    std::cout << "extension_20_time left alignment time: " << extension_20_time << "us" << std::endl;
//    std::cout << "extension_30_time left alignment time: " << extension_30_time << "us" << std::endl;
//    std::cout << "extension_100_time left alignment time: " << extension_100_time << "us" << std::endl;
//}
