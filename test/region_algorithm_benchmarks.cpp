//
//  region_algorithm_benchmarks.cpp
//  Octopus
//
//  Created by Daniel Cooke on 14/07/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

//#include <iostream>
//#include <algorithm>
//
//#include "benchmark_utils.hpp"
//#include "test_common.hpp"
//#include "genomic_region.hpp"
//#include "mappable_algorithms.hpp"
//#include "read_manager.hpp"

//BOOST_AUTO_TEST_CASE(overlap_range performance)
//{
//    ReadManager a_read_manager(std::vector<std::string> {human_1000g_bam1});
//    
//    auto sample = a_read_manager.get_sample_ids().front();
//    
//    GenomicRegion region {"1", 0, 249250621};
//    
//    auto reads = a_read_manager.fetch_reads(sample, region);
//    
//    GenomicRegion test_region {"1", 10000000, 15000000};
//    
//    std::size_t c1 {}, c2 {};
//    
//    std::cout << "num reads " << reads.size() << std::endl;
//    
//    auto max_read_size = size(*largest(reads.cbegin(), reads.cend()));
//    
//    auto f_overlap_range1 = [&reads, &test_region, &c1] () {
//        auto overlapped = overlap_range(reads.cbegin(), reads.cend(), test_region, true);
//        c1 = std::count_if(overlapped.begin(), overlapped.end(), [&test_region] (const auto& read) { return overlaps(read, test_region); });
//    };
//    
//    auto f_overlap_range2 = [&reads, &test_region, &c1, max_read_size] () {
//        auto overlapped = overlap_range(reads.cbegin(), reads.cend(), test_region, max_read_size);
//        c1 = std::count_if(overlapped.begin(), overlapped.end(), [&test_region] (const auto& read) { return overlaps(read, test_region); });
//    };
//    
//    auto f_overlap_range3 = [&reads, &test_region, &c1] () {
//        auto overlapped = overlap_range(reads.cbegin(), reads.cend(), test_region);
//        c1 = std::count_if(overlapped.begin(), overlapped.end(), [&test_region] (const auto& read) { return overlaps(read, test_region); });
//    };
//    
////    auto f_read_manager_fetch = [&a_read_manager, &sample, &test_region, &c2] () {
////        auto test_reads = a_read_manager.fetch_reads(sample, test_region);
////        c2 = std::count_if(test_reads.begin(), test_reads.end(), [&test_region] (const auto& read) { return overlaps(read, test_region); });
////    };
//    
//    auto time1 = benchmark<std::chrono::microseconds>(f_overlap_range1, 10).count();
//    auto time2 = benchmark<std::chrono::microseconds>(f_overlap_range2, 10).count();
//    auto time3 = benchmark<std::chrono::microseconds>(f_overlap_range3, 10).count();
//    //auto time2 = benchmark<std::chrono::microseconds>(f_read_manager_fetch, 10).count();
//    
//    std::cout << c1 << std::endl;
//    //std::cout << c2 << std::endl;
//    
//    std::cout << "overlap_range1 time " << time1 << std::endl;
//    std::cout << "overlap_range2 time " << time2 << std::endl;
//    std::cout << "overlap_range3 time " << time3 << std::endl;
//    //std::cout << "read_manager_fetch time " << time2 << std::endl;
//}
