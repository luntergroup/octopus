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

#include "test_common.hpp"
#include "genomic_region.hpp"
#include "aligned_read.hpp"
#include "htslib_sam_facade.hpp"
#include "read_manager.hpp"
#include "mock_objects.hpp"
#include "mappable_algorithms.hpp"

using std::cout;
using std::endl;

BOOST_AUTO_TEST_SUITE(Components)

BOOST_AUTO_TEST_CASE(aligned_read_copies_and_moves_correctly)
{
    AlignedRead read {get_mock_region(), "ACGT", AlignedRead::Qualities {1, 2, 3, 4},
        *parse_cigar_string("4M"), 10, AlignedRead::Flags {}, "1", 10, 30,
        AlignedRead::NextSegment::Flags {}};
    
    BOOST_CHECK(read.is_chimeric());
    BOOST_CHECK(read.get_next_segment().get_inferred_template_length() == 30);
    
    auto moved_read = std::move(read);
    
    BOOST_CHECK(moved_read.is_chimeric());
    BOOST_CHECK(moved_read.get_next_segment().get_inferred_template_length() == 30);
    
    auto copied_read = moved_read;
    
    BOOST_CHECK(copied_read.is_chimeric());
    BOOST_CHECK(copied_read.get_next_segment().get_inferred_template_length() == 30);
    BOOST_CHECK(moved_read == copied_read);
}

BOOST_AUTO_TEST_CASE(aligned_read_overlap_sanity_checks)
{
    BOOST_REQUIRE(test_file_exists(NA12878_low_coverage));
    
    ReadManager read_manager {NA12878_low_coverage};
    
    GenomicRegion region {"4", 93235280, 93235585};
    
    const auto sample = read_manager.get_samples().front();
    
    const auto reads = read_manager.fetch_reads(sample, region);
    
    BOOST_CHECK(reads.size() == 9);
    
    BOOST_REQUIRE(std::is_sorted(std::cbegin(reads), std::cend(reads)));
    
    GenomicRegion sub_region1 {"4", 93235300, 93235345};
    GenomicRegion sub_region2 {"4", 93235390, 93235445};
    GenomicRegion sub_region3 {"4", 93235445, 93235490};
    GenomicRegion sub_region4 {"4", 93235475, 93235490};
    GenomicRegion sub_region5 {"4", 93235500, 93235515};
    GenomicRegion sub_region6 {"4", 93235575, 93235580};
    
    auto reads_in_sub_region1 = overlap_range(reads, sub_region1);
    auto reads_in_sub_region2 = overlap_range(reads, sub_region2);
    auto reads_in_sub_region3 = overlap_range(reads, sub_region3);
    auto reads_in_sub_region4 = overlap_range(reads, sub_region4);
    auto reads_in_sub_region5 = overlap_range(reads, sub_region5);
    auto reads_in_sub_region6 = overlap_range(reads, sub_region6);
    
    BOOST_CHECK(size(reads_in_sub_region1) == 3);
    BOOST_CHECK(size(reads_in_sub_region2) == 2);
    BOOST_CHECK(size(reads_in_sub_region3) == 5);
    BOOST_CHECK(size(reads_in_sub_region4) == 4);
    BOOST_CHECK(size(reads_in_sub_region5) == 4);
    BOOST_CHECK(size(reads_in_sub_region6) == 1);
    
    BOOST_CHECK(count_shared(reads, sub_region1, sub_region2) == 0);
    BOOST_CHECK(count_shared(reads, sub_region2, sub_region3) == 2);
    BOOST_CHECK(count_shared(reads, sub_region2, sub_region4) == 1);
    BOOST_CHECK(count_shared(reads, sub_region2, sub_region5) == 0);
    BOOST_CHECK(count_shared(reads, sub_region3, sub_region4) == 4);
    BOOST_CHECK(count_shared(reads, sub_region2, sub_region5) == 0);
    BOOST_CHECK(count_shared(reads, sub_region3, sub_region5) == 3);
    BOOST_CHECK(count_shared(reads, sub_region4, sub_region5) == 3);
    BOOST_CHECK(count_shared(reads, sub_region4, sub_region6) == 0);
    BOOST_CHECK(count_shared(reads, sub_region5, sub_region6) == 1);
    
    BOOST_CHECK(!has_shared(reads, sub_region1, sub_region2));
    BOOST_CHECK(has_shared(reads, sub_region2, sub_region3));
    BOOST_CHECK(has_shared(reads, sub_region2, sub_region4));
    BOOST_CHECK(!has_shared(reads, sub_region2, sub_region5));
    BOOST_CHECK(has_shared(reads, sub_region3, sub_region4));
    BOOST_CHECK(!has_shared(reads, sub_region2, sub_region5));
    BOOST_CHECK(has_shared(reads, sub_region3, sub_region5));
    BOOST_CHECK(has_shared(reads, sub_region4, sub_region5));
    BOOST_CHECK(!has_shared(reads, sub_region4, sub_region6));
    BOOST_CHECK(has_shared(reads, sub_region5, sub_region6));
}

BOOST_AUTO_TEST_CASE(can_splice_CigarString)
{
    const auto cigar = *parse_cigar_string("5M1D10M3I4M");
    
    BOOST_CHECK(splice(cigar, 3, 10)  == *parse_cigar_string("2M1D7M"));
    BOOST_CHECK(splice(cigar, 3, 15)  == *parse_cigar_string("2M1D10M2I"));
    BOOST_CHECK(splice(cigar, 0, 10)  == *parse_cigar_string("5M1D4M"));
    BOOST_CHECK(splice(cigar, 0, 50)  == cigar);
    BOOST_CHECK(splice(cigar, 20, 10) == *parse_cigar_string("3M"));
    BOOST_CHECK(splice(cigar, 20, 3)  == *parse_cigar_string("3M"));
    BOOST_CHECK(splice(cigar, 24, 10) == *parse_cigar_string(""));
    BOOST_CHECK(splice(cigar, 16, 7)  == *parse_cigar_string("3I4M"));
}

BOOST_AUTO_TEST_CASE(can_splice_reads)
{
    AlignedRead read {
        GenomicRegion {"1", 100, 120},
        "AAAAACCCCCCCCCCGGGTTTT",
        AlignedRead::Qualities(23, 0),
        *parse_cigar_string("5M1D10M3I4M"),
        0,
        AlignedRead::Flags {}
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
    BOOST_REQUIRE(test_file_exists(NA12878_low_coverage));
    
    ReadManager a_read_manager {NA12878_low_coverage};
    
    GenomicRegion a_region {"4", 93235280, 93235585};
    
    auto sample_ids = a_read_manager.get_samples();
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
