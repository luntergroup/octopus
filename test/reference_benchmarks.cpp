//
//  reference_benchmarks.cpp
//  Octopus
//
//  Created by Daniel Cooke on 11/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include <iostream>
#include <string>
#include <chrono>
#include <fstream>
#include "benchmark_utils.h"

#include "test_common.h"
#include "reference_genome.h"
#include "bioio.h"

//BOOST_AUTO_TEST_CASE(reference_benchmark)
//{
//    std::ifstream fasta(human_reference_fasta);
//    auto index = bioio::read_fasta_index(human_reference_fasta_index)["5"];
//    auto f_bioio = [&fasta, &index] () {
//        bioio::read_fasta_contig(fasta, index, 100000, 100);
//    };
//    
//    auto without_vptr = benchmark<std::chrono::nanoseconds>(f_bioio, 10).count();
//    
//    fasta.close();
//    
//    ReferenceGenomeFactory a_factory {};
//    ReferenceGenome reference(a_factory.make(human_reference_fasta));
//    
//    auto f_ref = [&reference] () {
//        reference.get_sequence(GenomicRegion("5", 100000, 100100));
//    };
//    
//    auto with_vptr = benchmark<std::chrono::nanoseconds>(f_ref, 10).count();
//    
//    std::cout << "Without vptr: " << without_vptr << "ns" << std::endl;
//    std::cout << "With vptr: " << with_vptr << "ns" << std::endl;
//    
//    BOOST_CHECK(true);
//}

//BOOST_AUTO_TEST_CASE(caching_reference_benchmark)
//{
//    ReferenceGenomeFactory a_factory {};
//    
//    ReferenceGenome human_normal {a_factory.make(human_reference_fasta)};
//    ReferenceGenome human_cached {a_factory.make(human_reference_fasta, 1000000)};
//    
//    GenomicRegion region1 {"1", 1000000, 1000100};
//    GenomicRegion region2 {"1", 1000030, 1000080};
//    GenomicRegion region3 {"1",  999990, 1000110};
//    GenomicRegion region4 {"1",  1000200, 1000300};
//    GenomicRegion region5 {"1",  1000300, 1000400};
//    GenomicRegion region6 {"1",  1000350, 1000450};
//    
//    std::vector<GenomicRegion> regions {region1, region2, region3, region4, region5, region6};
//    
//    human_cached.get_sequence(region1);
//    human_cached.get_sequence(region2);
//    human_cached.get_sequence(region3);
//    human_cached.get_sequence(region4);
//    human_cached.get_sequence(region5);
//    human_cached.get_sequence(region6);
//    
//    auto f_normal = [&human_normal, &regions] () {
//        for (const auto& region : regions) {
//            human_normal.get_sequence(region);
//        }
//    };
//    
//    auto f_cached = [&human_cached, &regions] () {
//        for (const auto& region : regions) {
//            human_cached.get_sequence(region);
//        }
//    };
//    
//    auto normal = benchmark<std::chrono::nanoseconds>(f_normal, 10).count();
//    auto cached = benchmark<std::chrono::nanoseconds>(f_cached, 10).count();
//    
//    std::cout << "normal: " << normal << "ns" << std::endl;
//    std::cout << "cached: " << cached << "ns" << std::endl;
//}
