//
//  genotype_benchmarks.cpp
//  Octopus
//
//  Created by Daniel Cooke on 23/04/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

//#include <iostream>
//#include <string>
//
//#include "test_common.hpp"
//#include "benchmark_utils.hpp"
//#include "reference_genome.hpp"
//#include "read_manager.hpp"
//#include "variant.hpp"
//#include "candidate_variant_generator.hpp"
//#include "alignment_candidate_variant_generator.hpp"
//#include "haplotype.hpp"
//#include "genotype.hpp"
//#include "allele.hpp"
//#include "haplotype_tree.hpp"

//BOOST_AUTO_TEST_CASE(genotype hashing benchmark)
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
