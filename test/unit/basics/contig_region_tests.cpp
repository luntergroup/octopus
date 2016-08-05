// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include <boost/test/unit_test.hpp>

#include <iostream>

#include <basics/contig_region.hpp>

namespace octopus { namespace test {

BOOST_AUTO_TEST_SUITE(basics)
BOOST_AUTO_TEST_SUITE(contig_region)

BOOST_AUTO_TEST_CASE(overlaps_is_consistent)
{
    ContigRegion r1 {0, 0}, r2 {0, 1}, r3 {1, 1}, r4 {0, 2}, r5 {2, 2};
    
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

BOOST_AUTO_TEST_CASE(is_before_is_consistent)
{
    ContigRegion r1 {0, 0}, r2 {0, 1}, r3 {1, 1}, r4 {0, 2}, r5 {2, 2};
    
    BOOST_CHECK(is_before(r1, r2));
    BOOST_CHECK(!is_before(r2, r1));
    
    // more..
}

BOOST_AUTO_TEST_CASE(is_after_is_consistent)
{
    ContigRegion r1 {0, 0}, r2 {0, 1}, r3 {1, 1}, r4 {0, 2}, r5 {2, 2};
    
    BOOST_CHECK(is_after(r2, r1));
    BOOST_CHECK(!is_after(r1, r2));
    
    // more..
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

} // namespace test
} // namespace octopus
