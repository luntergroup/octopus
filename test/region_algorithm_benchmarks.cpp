//
//  region_algorithm_benchmarks.cpp
//  Octopus
//
//  Created by Daniel Cooke on 14/07/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "catch.hpp"

#include <iostream>
#include <algorithm>

#include "benchmark_utils.h"
#include "test_common.h"
#include "genomic_region.h"
#include "region_algorithms.h"
#include "read_manager.h"

//TEST_CASE("overlap_range performance", "[region_algorithms, benchmark]")
//{
//    ReadManager a_read_manager(std::vector<std::string> {human_1000g_bam1});
//    
//    auto sample = a_read_manager.get_sample_ids().front();
//    
//    GenomicRegion region {"1", 0, 20000000};
//    
//    auto reads = a_read_manager.fetch_reads(sample, region);
//    
//    GenomicRegion test_region {"1", 10000000, 15000000};
//    
//    std::size_t c1 {}, c2 {};
//    
//    auto f_overlap_range = [&reads, &test_region, &c1] () {
//        auto overlapped = overlap_range(reads.cbegin(), reads.cend(), test_region);
//        //auto m = std::max_element(overlapped.first, overlapped.second, [] (const auto& lhs, const auto& rhs) { return size(lhs) < size(rhs); });
//        c1 = std::count_if(overlapped.begin(), overlapped.end(), [&test_region] (const auto& read) { return overlaps(read, test_region); });
//        //std::cout << c << std::endl;
//        //std::cout << size(*m) << std::endl;
//        //auto num_reads = std::distance(overlapped.begin(), overlapped.end());
//    };
//    
//    auto time = benchmark<std::chrono::microseconds>(f_overlap_range, 10).count();
//    
//    std::cout << c1 << std::endl;
//    
//    std::cout << "overlap_range time " << time << std::endl;
//}
