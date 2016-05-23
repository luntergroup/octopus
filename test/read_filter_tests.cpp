//
//  read_filter_tests.cpp
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

//#include "test_common.hpp"
//#include "genomic_region.hpp"
//#include "read_manager.hpp"
//#include "read_filter.hpp"
//#include "read_filters.hpp"

BOOST_AUTO_TEST_SUITE(Components)

BOOST_AUTO_TEST_CASE(read_filter_test)
{
//    BOOST_REQUIRE(test_file_exists(NA12878_low_coverage));
//    
//    ReadManager read_manager {NA12878_low_coverage};
//    
//    auto sample_ids = read_manager.samples();
//    
//    auto sample = sample_ids.front();
//    
//    GenomicRegion region {"X", 1000000, 1010000};
//    
//    auto reads = read_manager.fetch_reads(sample, region);
//    
//    BOOST_REQUIRE(std::is_sorted(std::cbegin(reads), std::cend(reads)));
//    
//    BOOST_CHECK(reads.size() == 485);
//    
//    Octopus::ReadFilter<std::vector<AlignedRead>::iterator> read_filter {};
    
//    read_filter.register_filter(Octopus::ReadFilters::IsNotSecondaryAlignment());
//    read_filter.register_filter(Octopus::ReadFilters::IsGoodMappingQuality(20));
//    read_filter.register_filter(Octopus::ReadFilters::HasSufficientGoodQualityBases(20, 10));
//    read_filter.register_filter(Octopus::ReadFilters::FilterDuplicates());
    
//    std::vector<AlignedRead> good_reads {}, bad_reads {};
//    good_reads.reserve(reads.size());
//    bad_reads.reserve(reads.size());
//    
//    read_filter.partition_copy(std::make_move_iterator(reads.begin()), std::make_move_iterator(reads.end()),
//                               context_back_inserter(good_reads), context_back_inserter(bad_reads));
//    
//    // TODO: check these numbers are actually correct!
//    BOOST_CHECK(good_reads.size() == 436);
//    BOOST_CHECK(bad_reads.size() == 49);
}

BOOST_AUTO_TEST_SUITE_END()
