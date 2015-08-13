//
//  read_reader_test.cpp
//  Octopus
//
//  Created by Daniel Cooke on 11/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include <iostream>
#include <string>
#include <iterator>
#include <vector>

#include "test_common.h"
#include "genomic_region.h"
#include "aligned_read.h"
#include "htslib_sam_facade.h"
#include "read_manager.h"
#include "mock_objects.h"
#include "mappable_algorithms.h"

using std::cout;
using std::endl;

BOOST_AUTO_TEST_SUITE(Components)

BOOST_AUTO_TEST_CASE(read_reader_handles_BAM)
{
    HtslibSamFacade a_reader {human_1000g_bam1};
    
    GenomicRegion a_region {"10", 1000000, 1000100};
    GenomicRegion another_region {"3", 100000, 100100};
    
    auto reads      = a_reader.fetch_reads(a_region).begin()->second;
    auto more_reads = a_reader.fetch_reads(another_region).begin()->second;
    
    reads.insert(std::end(reads), std::make_move_iterator(std::begin(more_reads)),
                  std::make_move_iterator(std::end(more_reads)));
    
    BOOST_CHECK(reads.size() == 25);
}

BOOST_AUTO_TEST_CASE(read_reader_handles_CRAM)
{
    HtslibSamFacade a_reader {human_1000g_cram};
    
    GenomicRegion a_region {"10", 1000000, 1000100};
    GenomicRegion another_region {"3", 100000, 100100};
    
    auto reads      = a_reader.fetch_reads(a_region).begin()->second;
    
    auto more_reads = a_reader.fetch_reads(another_region).begin()->second;
    
    reads.insert(std::end(reads), std::make_move_iterator(std::begin(more_reads)),
                 std::make_move_iterator(std::end(more_reads)));
    
    BOOST_CHECK(reads.size() == 25);
}

BOOST_AUTO_TEST_CASE(ReadManager_works_with_single_file)
{
    ReadManager a_read_manager(std::vector<std::string> {human_1000g_bam1});
    
    BOOST_CHECK(a_read_manager.get_num_samples() == 1);
    
    auto sample_ids = a_read_manager.get_sample_ids();
    auto the_sample_id = sample_ids.front();
    
    GenomicRegion a_big_region {"1", 9990, 10000};
    GenomicRegion a_small_region {"10", 1000000, 1000100};
    GenomicRegion another_region {"3", 100000, 100100};
    
    auto some_reads = a_read_manager.fetch_reads(the_sample_id, a_big_region);
    auto reads      = a_read_manager.fetch_reads(the_sample_id, a_small_region);
    auto more_reads = a_read_manager.fetch_reads(the_sample_id, another_region);
    
    reads.insert(std::end(reads), std::make_move_iterator(std::begin(more_reads)),
                 std::make_move_iterator(std::end(more_reads)));
    
    BOOST_CHECK(reads.size() == 25);
}

//BOOST_AUTO_TEST_CASE(read_manager_multiple_files_below_max_file_limit_test)
//{
//    std::vector<std::string> read_paths {human_1000g_bam1, human_1000g_bam2};
//    
//    unsigned max_open_files {2};
//    
//    ReadManager a_read_manager(read_paths, max_open_files);
//    
//    BOOST_CHECK(a_read_manager.get_num_samples() == 2);
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
//    BOOST_CHECK(big_region_reads1.size() == 61225);
//    BOOST_CHECK(small_region_reads1.size() == 10);
//    BOOST_CHECK(big_region_reads2.size() == 142606);
//    BOOST_CHECK(small_region_reads2.size() == 29);
//}

//BOOST_AUTO_TEST_CASE(read_manager_multiple_files_above_max_file_limit_test)
//{
//    std::vector<std::string> read_paths {human_1000g_bam1, human_1000g_bam2, human_1000g_bam3};
//    
//    unsigned max_open_files {2};
//    
//    ReadManager a_read_manager(read_paths, max_open_files);
//    
//    BOOST_CHECK(a_read_manager.get_num_samples() == 3);
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
//    BOOST_CHECK(big_reads1.size() == 142606);
//    BOOST_CHECK(big_reads2.size() == 61225);
//    BOOST_CHECK(big_reads3.size() == 56433);
//    
//    auto small_reads1 = a_read_manager.fetch_reads(first_sample, a_small_region);
//    auto small_reads2 = a_read_manager.fetch_reads(second_sample, a_small_region);
//    auto small_reads3 = a_read_manager.fetch_reads(third_sample, a_small_region);
//    
//    BOOST_CHECK(small_reads1.size() == 29);
//    BOOST_CHECK(small_reads2.size() == 10);
//    BOOST_CHECK(small_reads3.size() == 8);
//}

//BOOST_AUTO_TEST_CASE(different_file_type_test", "[read_manager]")
//{
//    std::vector<std::string> read_paths {human_1000g_bam2, human_1000g_cram};
//    
//    ReadManager a_read_manager(read_paths);
//    
//    BOOST_CHECK(a_read_manager.get_num_samples() == 2);
//    
//    GenomicRegion a_big_region {"1", 2000000, 3000000};
//    
//}

//BOOST_AUTO_TEST_CASE(files_without_all_chromosomes_test", "[read_manager]")
//{
//    std::vector<std::string> read_paths {human_1000g_bam1, human_1000g_bam1_chrom_20};
//    
//    ReadManager a_read_manager(read_paths);
//    
//    BOOST_CHECK(a_read_manager.get_num_samples() == 1);
//    
//    GenomicRegion a_common_region {"20",
//    GenomicRegion a {"1", 2000000, 3000000};
//    
//}

BOOST_AUTO_TEST_CASE(aligned_read_copies_and_moves_correctly)
{
    AlignedRead a_read {get_mock_region(), "ACGT", AlignedRead::Qualities {1, 2, 3, 4},
        parse_cigar_string("4M"), 10, AlignedRead::FlagData {}, "1", 10, 30,
        AlignedRead::NextSegment::FlagData {}};
    
    BOOST_CHECK(a_read.is_chimeric());
    BOOST_CHECK(a_read.get_next_segment()->get_inferred_template_length() == 30);
    
    auto a_moved_read = std::move(a_read);
    
    BOOST_CHECK(a_moved_read.is_chimeric());
    BOOST_CHECK(a_moved_read.get_next_segment()->get_inferred_template_length() == 30);
    
    auto a_copied_read = a_moved_read;
    
    BOOST_CHECK(a_copied_read.is_chimeric());
    BOOST_CHECK(a_copied_read.get_next_segment()->get_inferred_template_length() == 30);
    BOOST_CHECK(a_moved_read == a_copied_read);
}

BOOST_AUTO_TEST_CASE(aligned_read_overlap_sanity_checks)
{
    ReadManager a_read_manager(std::vector<std::string> {human_1000g_bam1});
    
    GenomicRegion a_region {"4", 93235280, 93235585};
    
    auto sample_ids = a_read_manager.get_sample_ids();
    auto the_sample_id = sample_ids.front();
    
    auto reads = a_read_manager.fetch_reads(the_sample_id, a_region);
    
    BOOST_CHECK(reads.size() == 9);
    
    std::sort(reads.begin(), reads.end());
    
    GenomicRegion sub_region1 {"4", 93235300, 93235345};
    GenomicRegion sub_region2 {"4", 93235390, 93235445};
    GenomicRegion sub_region3 {"4", 93235445, 93235490};
    GenomicRegion sub_region4 {"4", 93235475, 93235490};
    GenomicRegion sub_region5 {"4", 93235500, 93235515};
    GenomicRegion sub_region6 {"4", 93235575, 93235580};
    
    auto reads_in_sub_region1 = overlap_range(reads.cbegin(), reads.cend(), sub_region1);
    auto reads_in_sub_region2 = overlap_range(reads.cbegin(), reads.cend(), sub_region2);
    auto reads_in_sub_region3 = overlap_range(reads.cbegin(), reads.cend(), sub_region3);
    auto reads_in_sub_region4 = overlap_range(reads.cbegin(), reads.cend(), sub_region4);
    auto reads_in_sub_region5 = overlap_range(reads.cbegin(), reads.cend(), sub_region5);
    auto reads_in_sub_region6 = overlap_range(reads.cbegin(), reads.cend(), sub_region6);
    
    BOOST_CHECK(size(reads_in_sub_region1) == 3);
    BOOST_CHECK(size(reads_in_sub_region2) == 2);
    BOOST_CHECK(size(reads_in_sub_region3) == 5);
    BOOST_CHECK(size(reads_in_sub_region4) == 4);
    BOOST_CHECK(size(reads_in_sub_region5) == 4);
    BOOST_CHECK(size(reads_in_sub_region6) == 1);
    
    BOOST_CHECK(count_shared(reads.cbegin(), reads.cend(), sub_region1, sub_region2) == 0);
    BOOST_CHECK(count_shared(reads.cbegin(), reads.cend(), sub_region2, sub_region3) == 2);
    BOOST_CHECK(count_shared(reads.cbegin(), reads.cend(), sub_region2, sub_region4) == 1);
    BOOST_CHECK(count_shared(reads.cbegin(), reads.cend(), sub_region2, sub_region5) == 0);
    BOOST_CHECK(count_shared(reads.cbegin(), reads.cend(), sub_region3, sub_region4) == 4);
    BOOST_CHECK(count_shared(reads.cbegin(), reads.cend(), sub_region2, sub_region5) == 0);
    BOOST_CHECK(count_shared(reads.cbegin(), reads.cend(), sub_region3, sub_region5) == 3);
    BOOST_CHECK(count_shared(reads.cbegin(), reads.cend(), sub_region4, sub_region5) == 3);
    BOOST_CHECK(count_shared(reads.cbegin(), reads.cend(), sub_region4, sub_region6) == 0);
    BOOST_CHECK(count_shared(reads.cbegin(), reads.cend(), sub_region5, sub_region6) == 1);
    
    BOOST_CHECK(!has_shared(reads.cbegin(), reads.cend(), sub_region1, sub_region2));
    BOOST_CHECK(has_shared(reads.cbegin(), reads.cend(), sub_region2, sub_region3));
    BOOST_CHECK(has_shared(reads.cbegin(), reads.cend(), sub_region2, sub_region4));
    BOOST_CHECK(!has_shared(reads.cbegin(), reads.cend(), sub_region2, sub_region5));
    BOOST_CHECK(has_shared(reads.cbegin(), reads.cend(), sub_region3, sub_region4));
    BOOST_CHECK(!has_shared(reads.cbegin(), reads.cend(), sub_region2, sub_region5));
    BOOST_CHECK(has_shared(reads.cbegin(), reads.cend(), sub_region3, sub_region5));
    BOOST_CHECK(has_shared(reads.cbegin(), reads.cend(), sub_region4, sub_region5));
    BOOST_CHECK(!has_shared(reads.cbegin(), reads.cend(), sub_region4, sub_region6));
    BOOST_CHECK(has_shared(reads.cbegin(), reads.cend(), sub_region5, sub_region6));
}

BOOST_AUTO_TEST_CASE(can_splice_CigarString)
{
    auto cigar = parse_cigar_string("5M1D10M3I4M");
    
    BOOST_CHECK(operation_splice(cigar, 3, 10)  == parse_cigar_string("2M1D7M"));
    BOOST_CHECK(operation_splice(cigar, 3, 15)  == parse_cigar_string("2M1D10M2I"));
    BOOST_CHECK(operation_splice(cigar, 0, 10)  == parse_cigar_string("5M1D4M"));
    BOOST_CHECK(operation_splice(cigar, 0, 50)  == cigar);
    BOOST_CHECK(operation_splice(cigar, 20, 10) == parse_cigar_string("3M"));
    BOOST_CHECK(operation_splice(cigar, 20, 3)  == parse_cigar_string("3M"));
    BOOST_CHECK(operation_splice(cigar, 24, 10) == parse_cigar_string(""));
    BOOST_CHECK(operation_splice(cigar, 16, 7)  == parse_cigar_string("3I4M"));
}

BOOST_AUTO_TEST_CASE(can_splice_reads)
{
    AlignedRead read {
        GenomicRegion {"1", 100, 120},
        "AAAAACCCCCCCCCCGGGTTTT",
        AlignedRead::Qualities(23, 0),
        parse_cigar_string("5M1D10M3I4M"),
        0,
        AlignedRead::FlagData {}
    };
    
    BOOST_CHECK(splice(read, GenomicRegion {"1", 100, 105}).get_sequence() == "AAAAA");
    BOOST_CHECK(splice(read, GenomicRegion {"1", 100, 106}).get_sequence() == "AAAAA");
    BOOST_CHECK(splice(read, GenomicRegion {"1", 100, 107}).get_sequence() == "AAAAAC");
    BOOST_CHECK(splice(read, GenomicRegion {"1", 100, 110}).get_sequence() == "AAAAACCCC");
    BOOST_CHECK(splice(read, GenomicRegion {"1", 100, 116}).get_sequence() == "AAAAACCCCCCCCCC");
    BOOST_CHECK(splice(read, GenomicRegion {"1", 100, 117}).get_sequence() == "AAAAACCCCCCCCCCGGGT");
    BOOST_CHECK(splice(read, GenomicRegion {"1", 100, 118}).get_sequence() == "AAAAACCCCCCCCCCGGGTT");
    BOOST_CHECK(splice(read, GenomicRegion {"1", 100, 119}).get_sequence() == "AAAAACCCCCCCCCCGGGTTT");
    BOOST_CHECK(splice(read, GenomicRegion {"1", 100, 120}) == read);
}

BOOST_AUTO_TEST_CASE(AlignedRead_can_be_compressed_and_decompressed)
{
    ReadManager a_read_manager(std::vector<std::string> {human_1000g_bam1});
    
    GenomicRegion a_region {"4", 93235280, 93235585};
    
    auto sample_ids = a_read_manager.get_sample_ids();
    auto the_sample_id = sample_ids.front();
    
    auto reads = a_read_manager.fetch_reads(the_sample_id, a_region);
    
    auto reads_copy = reads;
    
    for (auto& read : reads_copy) {
        read.compress();
    }
    
    for (auto& read : reads_copy) {
        read.decompress();
    }
    
    BOOST_CHECK(reads == reads_copy);
}

BOOST_AUTO_TEST_SUITE_END()
