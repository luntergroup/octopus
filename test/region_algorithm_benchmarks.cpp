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

TEST_CASE("overlap_range performance", "[region_algorithms, benchmark]")
{
    ReadManager a_read_manager(std::vector<std::string> {human_1000g_bam1});
    
    auto sample = a_read_manager.get_sample_ids().front();
    
    GenomicRegion region {"1", 0, 20000000};
    
    auto reads = a_read_manager.fetch_reads(sample, region);
    
    GenomicRegion test_region {"1", 10000000, 15000000};
    
    auto f_overlap_range = [&reads, &test_region] () {
        auto overlapped = overlap_range(reads.cbegin(), reads.cend(), test_region);
        auto m = std::max_element(overlapped.first, overlapped.second, [] (const auto& lhs, const auto& rhs) { return size(lhs) < size(rhs); });
        //std::cout << size(*m) << std::endl;
        //auto num_reads = std::distance(overlapped.begin(), overlapped.end());
    };
    
    auto f_overlap_range2 = [&reads, &test_region] () {
        auto overlapped = overlap_range2(reads.cbegin(), reads.cend(), test_region);
        auto m = std::max_element(overlapped.begin(), overlapped.end(), [] (const auto& lhs, const auto& rhs) { return size(lhs) < size(rhs); });
        //std::cout << size(*m) << std::endl;
        //auto num_reads = std::distance(overlapped.begin(), overlapped.end());
    };
    
    auto f_overlap_range3 = [&reads, &test_region] () {
        auto overlapped = overlap_range3(reads.cbegin(), reads.cend(), test_region, true);
        auto m = std::max_element(overlapped.begin(), overlapped.end(), [] (const auto& lhs, const auto& rhs) { return size(lhs) < size(rhs); });
        //std::cout << size(*m) << std::endl;
        //auto num_reads = std::distance(overlapped.begin(), overlapped.end());
    };
    
    auto f_overlap_ranges = [&reads, &test_region] () {
        auto overlapped = overlap_ranges(reads.cbegin(), reads.cend(), test_region);
        std::vector<AlignedRead> maxes {};
        maxes.reserve(overlapped.size());
        std::transform(overlapped.cbegin(), overlapped.cend(), std::back_inserter(maxes), [] (const auto& range) { return *std::max_element(range.first, range.second, [] (const auto& lhs, const auto& rhs) { return size(lhs) < size(rhs); }); });
        auto m = std::max_element(maxes.begin(), maxes.end(), [] (const auto& lhs, const auto& rhs) { return size(lhs) < size(rhs); });
        //std::cout << size(*m) << std::endl;
    };
    
    auto time1 = benchmark<std::chrono::microseconds>(f_overlap_range, 10).count();
    auto time2 = benchmark<std::chrono::microseconds>(f_overlap_range2, 10).count();
    auto time3 = benchmark<std::chrono::microseconds>(f_overlap_range3, 10).count();
    auto time4 = benchmark<std::chrono::microseconds>(f_overlap_ranges, 10).count();
    
    std::cout << "overlap_range time " << time1 << std::endl;
    std::cout << "overlap_range2 time " << time2 << std::endl;
    std::cout << "overlap_range3 time " << time3 << std::endl;
    std::cout << "overlap_ranges time " << time4 << std::endl;
}
