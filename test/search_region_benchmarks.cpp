//
//  search_region_benchmarks.cpp
//  Octopus
//
//  Created by Daniel Cooke on 17/07/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include <iostream>
#include <string>
#include <iterator>
#include <vector>
#include <algorithm>
#include <set>
#include <unordered_map>

#include "test_common.h"
#include "mock_objects.h"
#include "benchmark_utils.h"
#include "reference_genome.h"
#include "mappable_algorithms.h"
#include "test_common.h"
#include "read_manager.h"
#include "read_filter.h"
#include "read_filters.h"
#include "read_utils.h"
#include "allele.h"
#include "variant.h"
#include "candidate_variant_generator.h"
#include "alignment_candidate_variant_generator.h"
#include "haplotype.h"
#include "genotype.h"
#include "haplotype_tree.h"
#include "read_model.h"
#include "search_regions.h"
#include "read_utils.h"

using std::cout;
using std::endl;

using Octopus::advance_region;
using Octopus::cover_region;

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
