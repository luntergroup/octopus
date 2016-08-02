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
#include "read_manager.hpp"
#include "mock_objects.hpp"
#include "mappable_algorithms.hpp"

using std::cout;
using std::endl;

namespace octopus { namespace test {

BOOST_AUTO_TEST_SUITE(Components)

BOOST_AUTO_TEST_CASE(aligned_read_copies_and_moves_correctly)
{
    AlignedRead read {get_mock_region(), "ACGT", AlignedRead::BaseQualityVector {1, 2, 3, 4},
        parse_cigar("4M"), 10, AlignedRead::Flags {}, "1", 10, 30,
        AlignedRead::Segment::Flags {}};
    
    BOOST_CHECK(read.has_other_segment());
    BOOST_CHECK(read.next_segment().inferred_template_length() == 30);
    
    const auto moved_read = std::move(read);
    
    BOOST_CHECK(moved_read.has_other_segment());
    BOOST_CHECK(moved_read.next_segment().inferred_template_length() == 30);
    
    const auto copied_read = moved_read;
    
    BOOST_CHECK(copied_read.has_other_segment());
    BOOST_CHECK(copied_read.next_segment().inferred_template_length() == 30);
    BOOST_CHECK(moved_read == copied_read);
}

BOOST_AUTO_TEST_CASE(aligned_read_works_with_mappable_algorithms)
{
    BOOST_REQUIRE(test_file_exists(NA12878_low_coverage));
    
    ReadManager read_manager {NA12878_low_coverage};
    
//    const GenomicRegion region {"4", 10'000'000, 10'100'000};
//    
//    const auto sample = read_manager.samples().front();
//    
//    const auto reads = read_manager.fetch_reads(sample, region);
//    
//    BOOST_REQUIRE(reads.size() == 5011);
//    
//    BOOST_REQUIRE(std::is_sorted(std::cbegin(reads), std::cend(reads)));
//    
//    const GenomicRegion sub_region1 {"4", 10'000'000, 10'000'100};
//    const GenomicRegion sub_region2 {"4", 10'001'000, 10'001'100};
//    const GenomicRegion sub_region3 {"4", 10'010'000, 10'010'100};
//    const GenomicRegion sub_region4 {"4", 10'050'000, 10'050'100};
//    const GenomicRegion sub_region5 {"4", 10'050'000, 10'050'500};
//    const GenomicRegion sub_region6 {"4", 10'099'900, 10'100'000};
//    
//    const auto reads_in_sub_region1 = overlap_range(reads, sub_region1);
//    const auto reads_in_sub_region2 = overlap_range(reads, sub_region2);
//    const auto reads_in_sub_region3 = overlap_range(reads, sub_region3);
//    const auto reads_in_sub_region4 = overlap_range(reads, sub_region4);
//    const auto reads_in_sub_region5 = overlap_range(reads, sub_region5);
//    const auto reads_in_sub_region6 = overlap_range(reads, sub_region6);
//    
//    BOOST_CHECK(size(reads_in_sub_region1) == 14);
//    BOOST_CHECK(size(reads_in_sub_region2) == 7);
//    BOOST_CHECK(size(reads_in_sub_region3) == 8);
//    BOOST_CHECK(size(reads_in_sub_region4) == 6);
//    BOOST_CHECK(size(reads_in_sub_region5) == 19);
//    BOOST_CHECK(size(reads_in_sub_region6) == 17);
//    
//    BOOST_CHECK(count_shared(reads, sub_region1, sub_region2) == 0);
//    BOOST_CHECK(count_shared(reads, sub_region2, sub_region3) == 2);
//    BOOST_CHECK(count_shared(reads, sub_region2, sub_region4) == 1);
//    BOOST_CHECK(count_shared(reads, sub_region2, sub_region5) == 0);
//    BOOST_CHECK(count_shared(reads, sub_region3, sub_region4) == 4);
//    BOOST_CHECK(count_shared(reads, sub_region2, sub_region5) == 0);
//    BOOST_CHECK(count_shared(reads, sub_region3, sub_region5) == 3);
//    BOOST_CHECK(count_shared(reads, sub_region4, sub_region5) == 3);
//    BOOST_CHECK(count_shared(reads, sub_region4, sub_region6) == 0);
//    BOOST_CHECK(count_shared(reads, sub_region5, sub_region6) == 1);
//    
//    BOOST_CHECK(!has_shared(reads, sub_region1, sub_region2));
//    BOOST_CHECK(has_shared(reads, sub_region2, sub_region3));
//    BOOST_CHECK(has_shared(reads, sub_region2, sub_region4));
//    BOOST_CHECK(!has_shared(reads, sub_region2, sub_region5));
//    BOOST_CHECK(has_shared(reads, sub_region3, sub_region4));
//    BOOST_CHECK(!has_shared(reads, sub_region2, sub_region5));
//    BOOST_CHECK(has_shared(reads, sub_region3, sub_region5));
//    BOOST_CHECK(has_shared(reads, sub_region4, sub_region5));
//    BOOST_CHECK(!has_shared(reads, sub_region4, sub_region6));
//    BOOST_CHECK(has_shared(reads, sub_region5, sub_region6));
}

BOOST_AUTO_TEST_CASE(can_splice_CigarString)
{
    const auto cigar = parse_cigar("5M1D10M3I4M");
    
    BOOST_CHECK(splice(cigar, 3, 10)  == parse_cigar("2M1D7M"));
    BOOST_CHECK(splice(cigar, 3, 15)  == parse_cigar("2M1D10M2I"));
    BOOST_CHECK(splice(cigar, 0, 10)  == parse_cigar("5M1D4M"));
    BOOST_CHECK(splice(cigar, 0, 50)  == cigar);
    BOOST_CHECK(splice(cigar, 20, 10) == parse_cigar("3M"));
    BOOST_CHECK(splice(cigar, 20, 3)  == parse_cigar("3M"));
    BOOST_CHECK(splice(cigar, 24, 10) == parse_cigar(""));
    BOOST_CHECK(splice(cigar, 16, 7)  == parse_cigar("3I4M"));
}

BOOST_AUTO_TEST_CASE(can_splice_reads)
{
    const AlignedRead read {
        GenomicRegion {"1", 100, 120},
        "AAAAACCCCCCCCCCGGGTTTT",
        AlignedRead::BaseQualityVector(23, 0),
        parse_cigar("5M1D10M3I4M"),
        0,
        AlignedRead::Flags {}
    };
    
    BOOST_CHECK(splice(read, GenomicRegion {"1", 100, 105}).sequence() == "AAAAA");
    BOOST_CHECK(splice(read, GenomicRegion {"1", 100, 106}).sequence() == "AAAAA");
    BOOST_CHECK(splice(read, GenomicRegion {"1", 100, 107}).sequence() == "AAAAAC");
    BOOST_CHECK(splice(read, GenomicRegion {"1", 100, 110}).sequence() == "AAAAACCCC");
    BOOST_CHECK(splice(read, GenomicRegion {"1", 100, 116}).sequence() == "AAAAACCCCCCCCCC");
    BOOST_CHECK(splice(read, GenomicRegion {"1", 100, 117}).sequence() == "AAAAACCCCCCCCCCGGGT");
    BOOST_CHECK(splice(read, GenomicRegion {"1", 100, 118}).sequence() == "AAAAACCCCCCCCCCGGGTT");
    BOOST_CHECK(splice(read, GenomicRegion {"1", 100, 119}).sequence() == "AAAAACCCCCCCCCCGGGTTT");
    BOOST_CHECK(splice(read, GenomicRegion {"1", 100, 120}) == read);
}

BOOST_AUTO_TEST_CASE(AlignedRead_can_be_compressed_and_decompressed)
{
    BOOST_REQUIRE(test_file_exists(NA12878_low_coverage));
    
    ReadManager read_manager {NA12878_low_coverage};
    
    const GenomicRegion region {"4", 93235280, 93235585};
    
    const auto sample = read_manager.samples().front();
    
    auto reads = read_manager.fetch_reads(sample, region);
    
    BOOST_CHECK(!reads.empty());
    
    const auto reads_copy = reads;
    
//    for (auto& read : reads) {
//        read.compress();
//    }
//    
//    for (auto& read : reads) {
//        read.decompress();
//    }
    
    BOOST_CHECK(reads == reads_copy);
}

BOOST_AUTO_TEST_SUITE_END()

} // namespace test
} // namespace octopus
