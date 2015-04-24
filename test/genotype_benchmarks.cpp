//
//  genotype_benchmarks.cpp
//  Octopus
//
//  Created by Daniel Cooke on 23/04/2015.
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
#include "candidate_variant_generator.h"
#include "alignment_candidate_variant_generator.h"
#include "haplotype.h"
#include "genotype.h"
#include "allele.h"
#include "haplotype_tree.h"

//TEST_CASE("genotype hashing benchmark", "[genotype,benchmark]")
//{
//    ReferenceGenomeFactory a_factory {};
//    ReferenceGenome human {a_factory.make(human_reference_fasta)};
//    
//    auto a_region = parse_region("1:10000000-10001000", human);
//    
//    Haplotype haplotype {human, a_region};
//    
//    Genotype genotype {haplotype, haplotype};
//    
//    auto f_genotype_hash = [&genotype] () {
//        std::hash<Genotype>()(genotype);
//    };
//    
//    auto hash_time = benchmark<std::chrono::nanoseconds>(f_genotype_hash, 10000).count();
//    
//    std::cout << "hash_time: " << hash_time << "ns" << std::endl;
//}
