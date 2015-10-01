//
//  read_filter_benchmarks.cpp
//  Octopus
//
//  Created by Daniel Cooke on 07/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include <iostream>
#include <vector>
#include <algorithm> // std::sort
#include <iterator>  // std::back_inserter

#include "benchmark_utils.hpp"
#include "test_common.hpp"
#include "genomic_region.hpp"
#include "read_manager.hpp"
#include "read_filter.hpp"
#include "read_filters.hpp"

//BOOST_AUTO_TEST_CASE(read_filter_benchmarks)
//{
//    ReadManager a_read_manager(std::vector<std::string> {human_1000g_bam});
//    
//    auto sample_ids = a_read_manager.get_sample_ids();
//    
//    auto the_sample_id = sample_ids.at(0);
//    
//    GenomicRegion a_region {"1", 0, 20000000}; // ~1.3m reads
//    
//    std::vector<AlignedRead> some_reads = a_read_manager.fetch_reads(the_sample_id, a_region);
//    
//    if (!std::is_sorted(some_reads.begin(), some_reads.end())) {
//        std::sort(some_reads.begin(), some_reads.end());
//    }
//    
//    using ReadIterator = std::vector<AlignedRead>::const_iterator;
//    
//    ReadFilter<ReadIterator> a_read_filter {};
//    
//    // context-free filters
//    a_read_filter.register_filter(is_not_secondary_alignment);
//    a_read_filter.register_filter([] (const AlignedRead& the_read) {
//        return is_good_mapping_quality(the_read, 20);
//    });
//    a_read_filter.register_filter([] (const AlignedRead& the_read) {
//        return has_sufficient_good_quality_bases(the_read, 20, 10);
//    });
//    
//    // context-based filters
//    a_read_filter.register_filter(is_not_duplicate<ReadIterator>);
//    
//    auto f_filter = [&a_read_filter, &some_reads] () {
//        std::vector<AlignedRead> good_reads {}, bad_reads {};
//        good_reads.reserve(some_reads.size() / 2);
//        bad_reads.reserve(some_reads.size() / 2);
//        
//        a_read_filter.filter_reads(std::make_move_iterator(some_reads.begin()),
//                                   std::make_move_iterator(some_reads.end()),
//                                   ContextBackInserter(good_reads),
//                                   ContextBackInserter(bad_reads));
//    };
//    
//    auto filtering_time = benchmark<std::chrono::milliseconds>(f_filter, 1).count();
//    
//    std::cout << "filtering_time: " << filtering_time << "ms" << std::endl;
//}
