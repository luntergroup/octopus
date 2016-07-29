//
//  genome_walker_benchmarks.cpp
//  Octopus
//
//  Created by Daniel Cooke on 17/07/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

//#include <iostream>
//#include <string>
//#include <iterator>
//#include <vector>
//#include <algorithm>
//#include <set>
//#include <unordered_map>
//
//#include "test_common.hpp"
//#include "mock_objects.hpp"
//#include "benchmark_utils.hpp"
//#include "reference_genome.hpp"
//#include "mappable_algorithms.hpp"
//#include "test_common.hpp"
//#include "read_manager.hpp"
//#include "read_filter.hpp"
//#include "read_filters.hpp"
//#include "read_utils.hpp"
//#include "allele.hpp"
//#include "variant.hpp"
//#include "composer.hpp"
//#include "cigar_scanner.hpp"
//#include "haplotype.hpp"
//#include "genotype.hpp"
//#include "genome_walker.hpp"

//BOOST_AUTO_TEST_CASE(advance_region is reasonably fast)
//{
//    GenomicRegion::SizeType contig_size {10000000};
//    GenomicRegion::SizeType mean_read_size {100};
//    GenomicRegion::SizeType mean_variant_size {3};
//    std::size_t num_reads {5000000};
//    std::size_t num_variants {500000};
//    
//    auto candidates = generate_random_regions(contig_size, mean_variant_size, num_variants);
//    
//    GenomicRegion region {"test", 0, contig_size};
//    unsigned max_variants {10};
//    unsigned max_indicators {5};
//    
//    std::unordered_map<std::string, std::vector<GenomicRegion>> reads {};
//    reads.emplace("t", generate_random_regions(contig_size, mean_read_size, num_reads));
//    
//    cout << "starting cover benchmark" << endl;
//    
//    auto f_cover = [&candidates, &reads, &region, max_variants, max_indicators] () {
//        auto covers = cover_region(region, reads, candidates, max_variants, max_indicators);
//        cout << covers.size() << endl;
//    };
//    
//    auto time  = benchmark<std::chrono::microseconds>(f_cover, 10).count();
//    
//    cout << "cover time: " << time << "ms" << endl;
//}
