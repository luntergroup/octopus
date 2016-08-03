//
//  mappable_tests.cpp
//  Octopus
//
//  Created by Daniel Cooke on 02/02/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include <iostream>
#include <vector>
#include <algorithm>
#include <numeric>

#include <concepts/mappable.hpp>
#include <basics/contig_region.hpp>

using std::cout;
using std::endl;

namespace octopus {

BOOST_AUTO_TEST_SUITE(Components)
BOOST_AUTO_TEST_SUITE(MappableTests)

BOOST_AUTO_TEST_CASE(overlaps_is_consistent)
{
    ContigRegion r1 {0, 0};
    ContigRegion r2 {0, 1};
    ContigRegion r3 {1, 1};
    ContigRegion r4 {0, 2};
    ContigRegion r5 {2, 2};
    
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

BOOST_AUTO_TEST_CASE(is_after_and_is_before_are_consistent)
{
    ContigRegion r1 {0, 0};
    ContigRegion r2 {0, 1};
    ContigRegion r3 {1, 1};
    ContigRegion r4 {0, 2};
    ContigRegion r5 {2, 2};
    
    BOOST_CHECK(is_before(r1, r2));
    BOOST_CHECK(!is_before(r2, r1));
    BOOST_CHECK(is_after(r2, r1));
    BOOST_CHECK(!is_after(r1, r2));
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

} // namespace octopus
