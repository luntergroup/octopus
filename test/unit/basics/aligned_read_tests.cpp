// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include <utility>

#include "basics/genomic_region.hpp"
#include "basics/cigar_string.hpp"
#include "basics/aligned_read.hpp"

namespace octopus { namespace test {

BOOST_AUTO_TEST_SUITE(basics)
BOOST_AUTO_TEST_SUITE(aligned_read)

BOOST_AUTO_TEST_CASE(equal_operator_compares_)
{
    
}

BOOST_AUTO_TEST_CASE(comparison_operators_satisfy_mappable_requirments)
{
    
}

BOOST_AUTO_TEST_CASE(can_be_default_constructed)
{
    BOOST_CHECK_NO_THROW(AlignedRead {});
}

AlignedRead make_mock_read()
{
    return AlignedRead {
            GenomicRegion {"1", 0, 4}, "ACGT", AlignedRead::BaseQualityVector {1, 2, 3, 4},
            parse_cigar("4M"), 10, AlignedRead::Flags {}, "1", 10, 30, AlignedRead::Segment::Flags {}
    };
}

BOOST_AUTO_TEST_CASE(can_be_copied)
{
    const auto read1 = make_mock_read();
    AlignedRead read2 {};
    BOOST_REQUIRE_NO_THROW(read2 = read1);
    BOOST_CHECK_EQUAL(read1, read2);
}

BOOST_AUTO_TEST_CASE(can_be_moved)
{
    auto read1 = make_mock_read();
    AlignedRead read2 {};
    BOOST_REQUIRE_NO_THROW(read2 = std::move(read1));
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
    BOOST_CHECK_EQUAL(splice(read, GenomicRegion {"1", 100, 105}).sequence(), "AAAAA");
    BOOST_CHECK_EQUAL(splice(read, GenomicRegion {"1", 100, 106}).sequence(), "AAAAA");
    BOOST_CHECK_EQUAL(splice(read, GenomicRegion {"1", 100, 107}).sequence(), "AAAAAC");
    BOOST_CHECK_EQUAL(splice(read, GenomicRegion {"1", 100, 110}).sequence(), "AAAAACCCC");
    BOOST_CHECK_EQUAL(splice(read, GenomicRegion {"1", 100, 116}).sequence(), "AAAAACCCCCCCCCC");
    BOOST_CHECK_EQUAL(splice(read, GenomicRegion {"1", 100, 117}).sequence(), "AAAAACCCCCCCCCCGGGT");
    BOOST_CHECK_EQUAL(splice(read, GenomicRegion {"1", 100, 118}).sequence(), "AAAAACCCCCCCCCCGGGTT");
    BOOST_CHECK_EQUAL(splice(read, GenomicRegion {"1", 100, 119}).sequence(), "AAAAACCCCCCCCCCGGGTTT");
    BOOST_CHECK_EQUAL(splice(read, GenomicRegion {"1", 100, 120}), read);
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

} // namespace test
} // namespace octopus
