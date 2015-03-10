//
//  read_reader_test.cpp
//  Octopus
//
//  Created by Daniel Cooke on 11/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "catch.hpp"

#include <iostream>
#include <string>
#include <iterator>
#include <vector>

#include "test_common.h"
#include "genomic_region.h"
#include "htslib_read_facade.h"
#include "read_manager.h"
#include "mock_objects.h"

TEST_CASE("read_reader_open_test", "[read_reader]")
{
    HtslibReadFacade a_reader {human_1000g_bam};
    
    GenomicRegion a_region {"10", 1000000, 1000100};
    GenomicRegion another_region {"3", 100000, 100100};
    
    auto reads      = a_reader.fetch_reads(a_region).begin()->second;
    auto more_reads = a_reader.fetch_reads(another_region).begin()->second;
    
    reads.insert(std::end(reads), std::make_move_iterator(std::begin(more_reads)),
                  std::make_move_iterator(std::end(more_reads)));
    
    REQUIRE(reads.size() == 25);
}

TEST_CASE("read_manager_test", "[read_manager]")
{
    ReadManager a_read_manager(std::vector<std::string> {human_1000g_bam});
    
    auto sample_ids = a_read_manager.get_sample_ids();
    
    REQUIRE(sample_ids.size() == 1);
    
    auto the_sample_id = sample_ids.at(0);
    
    GenomicRegion a_big_region {"1", 9990, 10000};
    GenomicRegion a_small_region {"10", 1000000, 1000100};
    GenomicRegion another_region {"3", 100000, 100100};
    
    auto some_reads = a_read_manager.fetch_reads(the_sample_id, a_big_region);
    auto reads      = a_read_manager.fetch_reads(the_sample_id, a_small_region);
    auto more_reads = a_read_manager.fetch_reads(the_sample_id, another_region);
    
    reads.insert(std::end(reads), std::make_move_iterator(std::begin(more_reads)),
                 std::make_move_iterator(std::end(more_reads)));
    
    REQUIRE(reads.size() == 25);
}

TEST_CASE("read_copy_test", "[reads]")
{
    AlignedRead a_read {get_mock_region(), "ACGT", AlignedRead::Qualities {1, 2, 3, 4},
        parse_cigar_string("4M"), 10, "1", 10, 30, AlignedRead::SupplementaryData {},
        AlignedRead::MatePair::SupplementaryData {}};
    
    REQUIRE(a_read.has_mate_pair());
    REQUIRE(a_read.get_mate_pair()->get_insert_size() == 10);
    
    auto a_moved_read = std::move(a_read);
    
    REQUIRE(a_moved_read.has_mate_pair());
    REQUIRE(a_moved_read.get_mate_pair()->get_insert_size() == 10);
    
    auto a_copied_read = a_moved_read;
    
    REQUIRE(a_copied_read.has_mate_pair());
    REQUIRE(a_copied_read.get_mate_pair()->get_insert_size() == 10);
    REQUIRE(a_moved_read == a_copied_read);
}
