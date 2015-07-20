//
//  mappable_set_tests.cpp
//  Octopus
//
//  Created by Daniel Cooke on 19/07/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "catch.hpp"

#include <iostream>
#include <chrono>
#include <algorithm>

#include "test_common.h"
#include "genomic_region.h"
#include "read_manager.h"
#include "mappable_set.h"

using std::cout;
using std::endl;

TEST_CASE("MappableSet works like std::vector", "MappableSet")
{
    ReadManager a_read_manager(std::vector<std::string> {human_1000g_bam1});
    
    auto sample = a_read_manager.get_sample_ids().front();
    
    GenomicRegion a_region {"1", 0, 1000000};
    
    std::vector<AlignedRead> some_reads = a_read_manager.fetch_reads(sample, a_region);
    
    auto start = std::chrono::system_clock::now();
    
    MappableSet<AlignedRead> reads(some_reads.begin(), some_reads.end());
    
    auto end = std::chrono::system_clock::now();
    
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    
    cout << "duration MappableSet: " << duration.count() << endl;
    
    cout << reads.size() << endl;
    
    GenomicRegion overlap_region {"1", 500000, 600000};
    
    start = std::chrono::system_clock::now();
    
    auto overlapped = reads.overlap_range(overlap_region);
    
    auto has_overlap = reads.has_overlapped(overlap_region);
    
    auto has_contain = reads.has_contained(overlap_region);
    
    auto contained = reads.contained_range(overlap_region);
    
    end = std::chrono::system_clock::now();
    
    duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    
    cout << "duration MappableSet::overlap_range " << duration.count() << endl;
    
    cout << size(overlapped) << endl;
    cout << size(contained) << endl;
}
