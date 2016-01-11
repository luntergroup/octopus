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

#include "test_common.hpp"
#include "genomic_region.hpp"
#include "read_manager.hpp"
#include "read_filter.hpp"
#include "read_filters.hpp"
#include "context_iterators.hpp"
//#include "read_filter_factory.hpp"

using Octopus::context_back_inserter;

BOOST_AUTO_TEST_SUITE(Components)

BOOST_AUTO_TEST_CASE(read_filter_test)
{
    ReadManager read_manager {NA12878_low_coverage};
    
    auto sample_ids = read_manager.get_samples();
    
    auto sample = sample_ids.front();
    
    GenomicRegion region {"X", 1000000, 1010000};
    
    auto reads = read_manager.fetch_reads(sample, region);
    
    if (!std::is_sorted(reads.begin(), reads.end())) {
        std::sort(reads.begin(), reads.end());
    }
    
    BOOST_CHECK(reads.size() == 485);
    
    Octopus::ReadFilter<std::vector<AlignedRead>::const_iterator> read_filter {};
    
    read_filter.register_filter(Octopus::ReadFilters::is_not_secondary_alignment());
    read_filter.register_filter(Octopus::ReadFilters::is_good_mapping_quality(20));
    read_filter.register_filter(Octopus::ReadFilters::has_sufficient_good_quality_bases(20, 10));
    read_filter.register_filter(Octopus::ReadFilters::is_not_duplicate());
    
    std::vector<AlignedRead> good_reads {}, bad_reads {};
    good_reads.reserve(reads.size());
    bad_reads.reserve(reads.size());
    
    read_filter.filter_reads(std::make_move_iterator(reads.begin()), std::make_move_iterator(reads.end()),
                             context_back_inserter(good_reads), context_back_inserter(bad_reads));
    
    // TODO: check these numbers are actually correct!
    BOOST_CHECK(good_reads.size() == 436);
    BOOST_CHECK(bad_reads.size() == 49);
}

BOOST_AUTO_TEST_SUITE_END()
