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
#include "read_utils.h"
#include "read_filters.h"
#include "context_iterators.h"
#include "mappable_map.h"

using std::cout;
using std::endl;

TEST_CASE("MappableSet works like std::vector", "MappableSet")
{
    ReadManager a_read_manager(std::vector<std::string> {human_1000g_bam1});
    
    auto sample = a_read_manager.get_sample_ids().front();
    
    GenomicRegion a_region {"1", 0, 2000000};
    
    std::vector<AlignedRead> some_reads = a_read_manager.fetch_reads(sample, a_region);
    
    auto start = std::chrono::system_clock::now();
    
    MappableSet<AlignedRead> reads(std::make_move_iterator(some_reads.begin()), std::make_move_iterator(some_reads.end()));
    
    some_reads.clear();
    some_reads.shrink_to_fit();
    
    auto end = std::chrono::system_clock::now();
    
    auto duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    
    cout << "duration MappableSet: " << duration.count() << endl;
    
    cout << "total num reads " << reads.size() << endl;
    
    GenomicRegion overlap_region {"1", 500000, 600000};
    
    start = std::chrono::system_clock::now();
    
//    using ReadIterator = MappableSet<AlignedRead>::const_iterator;
//    
//    ReadFilter<ReadIterator> read_filter {};
//    
//    read_filter.register_filter(is_not_secondary_alignment);
//    read_filter.register_filter([] (const AlignedRead& the_read) {
//        return is_good_mapping_quality(the_read, 10);
//    });
//    read_filter.register_filter([] (const AlignedRead& the_read) {
//        return has_sufficient_good_quality_bases(the_read, 10, 5);
//    });
//    read_filter.register_filter(is_not_duplicate<ReadIterator>);
//    
//    MappableSet<AlignedRead> good_reads {}, bad_reads {};
//    good_reads.reserve(some_reads.size());
//    bad_reads.reserve(some_reads.size());
//    
//    read_filter.filter_reads(std::make_move_iterator(some_reads.begin()),
//                             std::make_move_iterator(some_reads.end()),
//                             Octopus::ContextInserter(good_reads),
//                             Octopus::ContextInserter(bad_reads));
    
//    cout << mean_coverage(reads, overlap_region) << endl;
//    
//    auto samples = downsample(reads, 50, 10);
//    
//    cout << mean_coverage(samples, overlap_region) << endl;
    
//    cout << good_reads.size() << endl;
//    cout << bad_reads.size() << endl;
//    
//    cout << good_reads.count_overlapped(overlap_region) << endl;
//    cout << bad_reads.count_overlapped(overlap_region) << endl;
    
//    GenomicRegion shared_region {"1", 600010, 600200};
//    
//    MappableMap<int, AlignedRead> read_map {};
//    read_map.emplace(1, reads);
//    
//    cout << *leftmost_overlapped(read_map, overlap_region) << endl;
    
    //auto filtered = filter_reads(std::move(reads), read_filter);
    
//    auto overlapped = reads.overlap_range(overlap_region);
//    
//    
//    cout << reads.count_overlapped(overlap_region) << endl;
//    
//    cout << "overlapped size " << size(overlapped) << endl;
//    
//    reads.erase(overlapped);
//    
//    cout << "new num reads " << reads.size() << endl;
//    
//    auto new_overlapped = reads.overlap_range(overlap_region);
//    
//    cout << "overlapped size " << size(new_overlapped) << endl;
    
//    auto has_overlap = reads.has_overlapped(overlap_region);
//    
//    auto has_contain = reads.has_contained(overlap_region);
//    
//    auto contained = reads.contained_range(overlap_region);
    
    end = std::chrono::system_clock::now();
    
    duration = std::chrono::duration_cast<std::chrono::microseconds>(end - start);
    
    cout << "duration MappableSet::overlap_range " << duration.count() << endl;
    
//    cout << size(overlapped) << endl;
//    cout << size(contained) << endl;
}
