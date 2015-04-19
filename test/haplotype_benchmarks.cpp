//
//  haplotype_benchmarks.cpp
//  Octopus
//
//  Created by Daniel Cooke on 13/04/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "catch.hpp"

#include <iostream>
#include <string>

#include "test_common.h"
#include "benchmark_utils.h"
#include "reference_genome.h"
#include "reference_genome_factory.h"
#include "read_manager.h"
#include "variant.h"
#include "variant_utils.h"
#include "variant_candidate_generator.h"
#include "alignment_candidate_variant_generator.h"
#include "haplotype.h"
#include "genotype.h"

//TEST_CASE("haplotype hashing benchmark", "[haplotype]")
//{
//    ReferenceGenomeFactory a_factory {};
//    ReferenceGenome human {a_factory.make(human_reference_fasta)};
//    
//    auto a_region = parse_region("1:10000000-10001000", human);
//    
//    Haplotype haplotype {human, a_region};
//    
//    auto f_haplotype_hash = [&haplotype] () {
//        std::hash<Haplotype>()(haplotype);
//    };
//    
//    auto hash_time = benchmark<std::chrono::nanoseconds>(f_haplotype_hash, 10000).count();
//    
//    std::cout << "hash_time: " << hash_time << "ns" << std::endl;
//}
