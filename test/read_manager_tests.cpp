//
//  read_manager_tests.cpp
//  Octopus
//
//  Created by Daniel Cooke on 01/02/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include <iostream>
#include <string>
#include <iterator>
#include <vector>

#include <boost/filesystem.hpp>

#include "test_common.hpp"
#include "genomic_region.hpp"
#include "aligned_read.hpp"
#include "htslib_sam_facade.hpp"
#include "read_manager.hpp"
#include "mock_objects.hpp"
#include "mappable_algorithms.hpp"

using std::cout;
using std::endl;

namespace fs = boost::filesystem;

BOOST_AUTO_TEST_SUITE(Components)

BOOST_AUTO_TEST_CASE(HtslibSamFacade_can_extract_reads_from_BAM_files)
{
    BOOST_REQUIRE(test_file_exists(NA12878_low_coverage));
    
    HtslibSamFacade reader {NA12878_low_coverage};
    
    GenomicRegion region1 {"10", 1000000, 1000100};
    GenomicRegion region2 {"3", 100000, 100100};
    
    auto reads1 = reader.fetch_reads(region1).begin()->second;
    auto reads2 = reader.fetch_reads(region2).begin()->second;
    
    reads1.insert(std::begin(reads2), std::end(reads2));
    
    BOOST_CHECK(reads1.size() == 25);
}

BOOST_AUTO_TEST_CASE(HtslibSamFacade_can_extract_reads_from_CRAM_files)
{
    BOOST_REQUIRE(test_file_exists(HG00101_cram));
    
    HtslibSamFacade reader {HG00101_cram};
    
    GenomicRegion region1 {"10", 1000000, 1000100};
    GenomicRegion region2 {"3", 100000, 100100};
    
    auto reads1 = reader.fetch_reads(region1).begin()->second;
    auto reads2 = reader.fetch_reads(region2).begin()->second;
    
    reads1.insert(std::begin(reads2), std::end(reads2));
    
    BOOST_CHECK(reads1.size() == 25);
}

BOOST_AUTO_TEST_CASE(ReadManager_can_extract_reads_from_BAM_files)
{
    BOOST_REQUIRE(test_file_exists(NA12878_low_coverage));
    
    ReadManager read_manager {NA12878_low_coverage};
    
    BOOST_CHECK(read_manager.num_samples() == 1);
    
    auto sample_ids = read_manager.get_samples();
    auto the_sample_id = sample_ids.front();
    
    GenomicRegion a_big_region {"1", 9990, 10000};
    GenomicRegion a_small_region {"10", 1000000, 1000100};
    GenomicRegion another_region {"3", 100000, 100100};
    
    auto some_reads = read_manager.fetch_reads(the_sample_id, a_big_region);
    auto reads      = read_manager.fetch_reads(the_sample_id, a_small_region);
    auto more_reads = read_manager.fetch_reads(the_sample_id, another_region);
    
    reads.insert(std::make_move_iterator(std::begin(more_reads)),
                 std::make_move_iterator(std::end(more_reads)));
    
    BOOST_CHECK(reads.size() == 25);
}

BOOST_AUTO_TEST_CASE(read_manager_multiple_files_below_max_file_limit_test)
{
    BOOST_REQUIRE(test_file_exists(NA12878_low_coverage));
    BOOST_REQUIRE(test_file_exists(HG00101));
    
    std::vector<fs::path> read_paths {NA12878_low_coverage, HG00101};
    
    constexpr unsigned max_open_files {2};
    
    ReadManager read_manager(read_paths, max_open_files);
    
    BOOST_CHECK(read_manager.num_samples() == 2);
    
    auto samples = read_manager.get_samples();
    
    GenomicRegion a_big_region {"1", 2000000, 3000000};
    GenomicRegion a_small_region {"10", 1000000, 1000100};
    
    auto first_sample = samples.at(0);
    
    auto big_region_reads1   = read_manager.fetch_reads(first_sample, a_big_region);
    auto small_region_reads1 = read_manager.fetch_reads(first_sample, a_small_region);
    
    auto second_sample = samples.at(1);
    
    auto big_region_reads2   = read_manager.fetch_reads(second_sample, a_big_region);
    auto small_region_reads2 = read_manager.fetch_reads(second_sample, a_small_region);
    
    BOOST_CHECK(big_region_reads1.size() == 61225);
    BOOST_CHECK(small_region_reads1.size() == 10);
    BOOST_CHECK(big_region_reads2.size() == 142606);
    BOOST_CHECK(small_region_reads2.size() == 29);
}

BOOST_AUTO_TEST_CASE(read_manager_multiple_files_above_max_file_limit_test)
{
    BOOST_REQUIRE(test_file_exists(NA12878_low_coverage));
    BOOST_REQUIRE(test_file_exists(HG00101));
    BOOST_REQUIRE(test_file_exists(HG00102));
    
    std::vector<fs::path> read_paths {NA12878_low_coverage, HG00101, HG00102};
    
    constexpr unsigned max_open_files {2};
    
    ReadManager read_manager(read_paths, max_open_files);
    
    BOOST_CHECK(read_manager.num_samples() == 3);
    
    auto samples = read_manager.get_samples();
    
    GenomicRegion a_big_region {"1", 2000000, 3000000};
    GenomicRegion a_small_region {"10", 1000000, 1000100};
    
    auto big_reads1 = read_manager.fetch_reads(samples[0], a_big_region);
    auto big_reads2 = read_manager.fetch_reads(samples[1], a_big_region);
    auto big_reads3 = read_manager.fetch_reads(samples[2], a_big_region);
    
    BOOST_CHECK(big_reads1.size() == 142606);
    BOOST_CHECK(big_reads2.size() == 61225);
    BOOST_CHECK(big_reads3.size() == 56433);
    
    auto small_reads1 = read_manager.fetch_reads(samples[0], a_small_region);
    auto small_reads2 = read_manager.fetch_reads(samples[1], a_small_region);
    auto small_reads3 = read_manager.fetch_reads(samples[2], a_small_region);
    
    BOOST_CHECK(small_reads1.size() == 29);
    BOOST_CHECK(small_reads2.size() == 10);
    BOOST_CHECK(small_reads3.size() == 8);
}

BOOST_AUTO_TEST_SUITE_END()
