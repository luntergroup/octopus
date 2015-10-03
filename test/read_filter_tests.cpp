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

using Octopus::ContextBackInserter;

BOOST_AUTO_TEST_SUITE(Components)

BOOST_AUTO_TEST_CASE(read_filter_test)
{
    ReadManager a_read_manager {human_1000g_bam1};
    
    auto sample_ids = a_read_manager.get_samples();
    
    auto the_sample_id = sample_ids.at(0);
    
    GenomicRegion a_region {"X", 1000000, 1010000};
    
    auto some_reads = a_read_manager.fetch_reads(the_sample_id, a_region);
    
    if (!std::is_sorted(some_reads.begin(), some_reads.end())) {
        std::sort(some_reads.begin(), some_reads.end());
    }
    
    BOOST_CHECK(some_reads.size() == 669);
    
    using ReadIterator = std::vector<AlignedRead>::const_iterator;
    
    Octopus::ReadFilter<ReadIterator> read_filter {};
    
    read_filter.register_filter(Octopus::ReadFilters::is_not_secondary_alignment());
    read_filter.register_filter(Octopus::ReadFilters::is_good_mapping_quality(20));
    read_filter.register_filter(Octopus::ReadFilters::has_sufficient_good_quality_bases(20, 10));
    read_filter.register_filter(Octopus::ReadFilters::is_not_duplicate());
    
    std::vector<AlignedRead> good_reads {}, bad_reads {};
    good_reads.reserve(some_reads.size());
    bad_reads.reserve(some_reads.size());
    
    read_filter.filter_reads(std::make_move_iterator(some_reads.begin()),
                             std::make_move_iterator(some_reads.end()),
                             ContextBackInserter(good_reads),
                             ContextBackInserter(bad_reads));
    
    // TODO: check these numbers are actually correct!
    BOOST_CHECK(good_reads.size() == 649);
    BOOST_CHECK(bad_reads.size() == 20);
}

BOOST_AUTO_TEST_SUITE_END()
