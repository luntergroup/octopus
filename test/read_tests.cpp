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
    HtslibReadFacade a_reader {human_1000g_bam1};
    
    GenomicRegion a_region {"10", 1000000, 1000100};
    GenomicRegion another_region {"3", 100000, 100100};
    
    auto reads      = a_reader.fetch_reads(a_region).begin()->second;
    auto more_reads = a_reader.fetch_reads(another_region).begin()->second;
    
    reads.insert(std::end(reads), std::make_move_iterator(std::begin(more_reads)),
                  std::make_move_iterator(std::end(more_reads)));
    
    REQUIRE(reads.size() == 25);
}

TEST_CASE("read_reader_open_test_cram", "[read_reader]")
{
    HtslibReadFacade a_reader {human_1000g_cram};
    
    GenomicRegion a_region {"10", 1000000, 1000100};
    GenomicRegion another_region {"3", 100000, 100100};
    
    auto reads      = a_reader.fetch_reads(a_region).begin()->second;
    
    auto more_reads = a_reader.fetch_reads(another_region).begin()->second;
    
    reads.insert(std::end(reads), std::make_move_iterator(std::begin(more_reads)),
                 std::make_move_iterator(std::end(more_reads)));
    
    REQUIRE(reads.size() == 25);
}

TEST_CASE("read_manager_single_file_test", "[read_manager]")
{
    ReadManager a_read_manager(std::vector<std::string> {human_1000g_bam1});
    
    REQUIRE(a_read_manager.get_num_samples() == 1);
    
    auto sample_ids = a_read_manager.get_sample_ids();
    
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

//TEST_CASE("read_manager_multiple_files_below_max_file_limit_test", "[read_manager")
//{
//    std::vector<std::string> read_paths {human_1000g_bam1, human_1000g_bam2};
//    
//    unsigned max_open_files {2};
//    
//    ReadManager a_read_manager(read_paths, max_open_files);
//    
//    REQUIRE(a_read_manager.get_num_samples() == 2);
//    
//    auto sample_ids = a_read_manager.get_sample_ids();
//    
//    GenomicRegion a_big_region {"1", 2000000, 3000000};
//    GenomicRegion a_small_region {"10", 1000000, 1000100};
//    
//    auto first_sample = sample_ids.at(0);
//    
//    auto big_region_reads1   = a_read_manager.fetch_reads(first_sample, a_big_region);
//    auto small_region_reads1 = a_read_manager.fetch_reads(first_sample, a_small_region);
//    
//    auto second_sample = sample_ids.at(1);
//    
//    auto big_region_reads2   = a_read_manager.fetch_reads(second_sample, a_big_region);
//    auto small_region_reads2 = a_read_manager.fetch_reads(second_sample, a_small_region);
//    
//    REQUIRE(big_region_reads1.size() == 61225);
//    REQUIRE(small_region_reads1.size() == 10);
//    REQUIRE(big_region_reads2.size() == 142606);
//    REQUIRE(small_region_reads2.size() == 29);
//}
//
//TEST_CASE("read_manager_multiple_files_above_max_file_limit_test", "[read_manager")
//{
//    std::vector<std::string> read_paths {human_1000g_bam1, human_1000g_bam2, human_1000g_bam3};
//    
//    unsigned max_open_files {2};
//    
//    ReadManager a_read_manager(read_paths, max_open_files);
//    
//    REQUIRE(a_read_manager.get_num_samples() == 3);
//    
//    auto sample_ids = a_read_manager.get_sample_ids();
//    
//    GenomicRegion a_big_region {"1", 2000000, 3000000};
//    GenomicRegion a_small_region {"10", 1000000, 1000100};
//    
//    auto first_sample = sample_ids.at(0);
//    auto second_sample = sample_ids.at(1);
//    auto third_sample = sample_ids.at(2);
//    
//    auto big_reads1 = a_read_manager.fetch_reads(first_sample, a_big_region);
//    auto big_reads2 = a_read_manager.fetch_reads(second_sample, a_big_region);
//    auto big_reads3 = a_read_manager.fetch_reads(third_sample, a_big_region);
//    
//    REQUIRE(big_reads1.size() == 142606);
//    REQUIRE(big_reads2.size() == 61225);
//    REQUIRE(big_reads3.size() == 56433);
//    
//    auto small_reads1 = a_read_manager.fetch_reads(first_sample, a_small_region);
//    auto small_reads2 = a_read_manager.fetch_reads(second_sample, a_small_region);
//    auto small_reads3 = a_read_manager.fetch_reads(third_sample, a_small_region);
//    
//    REQUIRE(small_reads1.size() == 29);
//    REQUIRE(small_reads2.size() == 10);
//    REQUIRE(small_reads3.size() == 8);
//}

//TEST_CASE("different_file_type_test", "[read_manager]")
//{
//    std::vector<std::string> read_paths {human_1000g_bam2, human_1000g_cram};
//    
//    ReadManager a_read_manager(read_paths);
//    
//    REQUIRE(a_read_manager.get_num_samples() == 2);
//    
//    GenomicRegion a_big_region {"1", 2000000, 3000000};
//    
//}

//TEST_CASE("files_without_all_chromosomes_test", "[read_manager]")
//{
//    std::vector<std::string> read_paths {human_1000g_bam1, human_1000g_bam1_chrom_20};
//    
//    ReadManager a_read_manager(read_paths);
//    
//    REQUIRE(a_read_manager.get_num_samples() == 1);
//    
//    GenomicRegion a_common_region {"20",
//    GenomicRegion a {"1", 2000000, 3000000};
//    
//}

TEST_CASE("read_copy_test", "[reads]")
{
    AlignedRead a_read {get_mock_region(), "ACGT", AlignedRead::Qualities {1, 2, 3, 4},
        parse_cigar_string("4M"), 10, AlignedRead::FlagData {}, "1", 10, 30,
        AlignedRead::NextSegment::FlagData {}};
    
    REQUIRE(a_read.is_chimeric());
    REQUIRE(a_read.get_next_segment()->get_inferred_template_length() == 30);
    
    auto a_moved_read = std::move(a_read);
    
    REQUIRE(a_moved_read.is_chimeric());
    REQUIRE(a_moved_read.get_next_segment()->get_inferred_template_length() == 30);
    
    auto a_copied_read = a_moved_read;
    
    REQUIRE(a_copied_read.is_chimeric());
    REQUIRE(a_copied_read.get_next_segment()->get_inferred_template_length() == 30);
    REQUIRE(a_moved_read == a_copied_read);
}
