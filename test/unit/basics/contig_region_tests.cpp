// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include <boost/test/unit_test.hpp>

#include <stdexcept>

#include <basics/contig_region.hpp>

namespace octopus { namespace test {

BOOST_AUTO_TEST_SUITE(basics)
BOOST_AUTO_TEST_SUITE(contig_region)

BOOST_AUTO_TEST_CASE(constructing_a_negative_region_is_an_error)
{
    BOOST_CHECK_NO_THROW((ContigRegion {0, 0}));
    BOOST_CHECK_NO_THROW((ContigRegion {0, 1}));
    BOOST_CHECK_THROW((ContigRegion {1, 0}), std::runtime_error);
}

BOOST_AUTO_TEST_CASE(ordering_is_by_begin_then_end)
{
    ContigRegion r1 {0, 0}, r2 {0, 1}, r3 {1, 1}, r4 {0, 2};
    
    BOOST_CHECK_NE(r1, r2);
    BOOST_CHECK_LT(r1, r2);
    
    BOOST_CHECK_NE(r2, r3);
    BOOST_CHECK_LT(r2, r3);
    
    BOOST_CHECK_NE(r1, r4);
    BOOST_CHECK_LT(r1, r4);
    
    BOOST_CHECK_NE(r2, r4);
    BOOST_CHECK_LT(r2, r4);
    
    BOOST_CHECK_NE(r3, r4);
    BOOST_CHECK_LT(r4, r3);
}

BOOST_AUTO_TEST_CASE(is_before_is_consistent)
{
    ContigRegion r1 {0, 0}, r2 {0, 1}, r3 {1, 1}, r4 {0, 2}, r5 {2, 2};
    
    BOOST_CHECK(!is_before(r1, r1));
    BOOST_CHECK(!is_before(r2, r2));
    
    BOOST_CHECK(is_before(r1, r2));
    BOOST_CHECK(!is_before(r2, r1));
    
    BOOST_CHECK(is_before(r1, r3));
    BOOST_CHECK(!is_before(r3, r1));
    
    BOOST_CHECK(is_before(r1, r4));
    BOOST_CHECK(!is_before(r4, r1));
    
    BOOST_CHECK(is_before(r4, r5));
    BOOST_CHECK(!is_before(r5, r4));
    
    BOOST_CHECK(!is_before(r3, r4));
    BOOST_CHECK(!is_before(r4, r3));
}

BOOST_AUTO_TEST_CASE(is_after_is_consistent)
{
    ContigRegion r1 {0, 0}, r2 {0, 1}, r3 {1, 1}, r4 {0, 2}, r5 {2, 2};
    
    BOOST_CHECK(!is_after(r1, r1));
    BOOST_CHECK(!is_after(r2, r2));
    
    BOOST_CHECK(is_after(r2, r1));
    BOOST_CHECK(!is_after(r1, r2));
    
    BOOST_CHECK(is_after(r3, r1));
    BOOST_CHECK(!is_after(r1, r3));
    
    BOOST_CHECK(is_after(r4, r1));
    BOOST_CHECK(!is_after(r1, r4));
    
    BOOST_CHECK(is_after(r5, r2));
    BOOST_CHECK(!is_after(r2, r5));
    
    BOOST_CHECK(is_after(r5, r3));
    BOOST_CHECK(!is_after(r3, r5));
    
    BOOST_CHECK(!is_after(r3, r4));
    BOOST_CHECK(!is_after(r3, r4));
}

BOOST_AUTO_TEST_CASE(overlap_size_returns_the_number_of_overlapped_positions)
{
    ContigRegion r1 {0, 0}, r2 {0, 1}, r3 {0, 2}, r4 {0, 4};
    
    BOOST_CHECK_EQUAL(overlap_size(r1, r1), 0);
    BOOST_CHECK_EQUAL(overlap_size(r1, r2), 0);
    BOOST_CHECK_EQUAL(overlap_size(r1, r3), 0);
    BOOST_CHECK_EQUAL(overlap_size(r1, r4), 0);
    BOOST_CHECK_EQUAL(overlap_size(r2, r1), 0);
    BOOST_CHECK_EQUAL(overlap_size(r3, r1), 0);
    BOOST_CHECK_EQUAL(overlap_size(r4, r1), 0);
    
    BOOST_CHECK_EQUAL(overlap_size(r2, r3), 1);
    BOOST_CHECK_EQUAL(overlap_size(r3, r2), 1);
    
    BOOST_CHECK_EQUAL(overlap_size(r2, r4), 1);
    BOOST_CHECK_EQUAL(overlap_size(r4, r2), 1);
    
    BOOST_CHECK_EQUAL(overlap_size(r3, r4), 2);
    BOOST_CHECK_EQUAL(overlap_size(r4, r3), 2);
}
    
BOOST_AUTO_TEST_CASE(overlaps_is_consistent)
{
    ContigRegion r1 {0, 0}, r2 {0, 1}, r3 {1, 1}, r4 {0, 2}, r5 {2, 2};
    
    BOOST_CHECK(overlaps(r1, r1)); // empty regions overlap
    BOOST_CHECK(overlaps(r2, r2));
    BOOST_CHECK(overlaps(r3, r3));
    BOOST_CHECK(overlaps(r4, r4));
    BOOST_CHECK(overlaps(r5, r5));
    
    BOOST_CHECK(overlaps(r1, r2));
    BOOST_CHECK(overlaps(r2, r1));
    BOOST_CHECK(!overlaps(r1, r3));
    BOOST_CHECK(!overlaps(r3, r1));
    BOOST_CHECK(overlaps(r2, r3));
    BOOST_CHECK(overlaps(r3, r2));
    
    BOOST_CHECK(overlaps(r1, r4));
    BOOST_CHECK(overlaps(r2, r4));
    BOOST_CHECK(overlaps(r3, r4));
    BOOST_CHECK(overlaps(r4, r1));
    BOOST_CHECK(overlaps(r4, r2));
    BOOST_CHECK(overlaps(r4, r3));
    
    BOOST_CHECK(!overlaps(r1, r5));
    BOOST_CHECK(!overlaps(r2, r5));
    BOOST_CHECK(!overlaps(r3, r5));
    BOOST_CHECK(!overlaps(r5, r1));
    BOOST_CHECK(!overlaps(r5, r2));
    BOOST_CHECK(!overlaps(r5, r3));
}

BOOST_AUTO_TEST_CASE(contains_is_consistent)
{
    ContigRegion r1 {0, 0}, r2 {0, 1}, r3 {1, 1}, r4 {0, 2}, r5 {2, 2};
    
    BOOST_CHECK(contains(r1, r1));
    BOOST_CHECK(contains(r2, r2));
    BOOST_CHECK(contains(r3, r3));
    BOOST_CHECK(contains(r4, r4));
    BOOST_CHECK(contains(r5, r5));
    
    BOOST_CHECK(contains(r2, r1));
    BOOST_CHECK(!contains(r1, r2));
    
    BOOST_CHECK(contains(r2, r3));
    BOOST_CHECK(!contains(r3, r2));
    
    BOOST_CHECK(contains(r4, r1));
    BOOST_CHECK(contains(r4, r2));
    BOOST_CHECK(contains(r4, r3));
    BOOST_CHECK(contains(r4, r5));
    BOOST_CHECK(!contains(r1, r4));
    BOOST_CHECK(!contains(r2, r4));
    BOOST_CHECK(!contains(r3, r4));
    BOOST_CHECK(!contains(r5, r4));
}

BOOST_AUTO_TEST_CASE(overlapping_empty_regions_are_considered_adjacent)
{
    ContigRegion r1 {0, 0}, r2 {0, 1}, r3 {1, 1}, r4 {0, 2}, r5 {2, 2};
    
    BOOST_CHECK(are_adjacent(r1, r1));
    BOOST_CHECK(are_adjacent(r3, r3));
    BOOST_CHECK(are_adjacent(r5, r5));
    
    BOOST_CHECK(!are_adjacent(r2, r2));
    BOOST_CHECK(!are_adjacent(r4, r4));
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

} // namespace test
} // namespace octopus
