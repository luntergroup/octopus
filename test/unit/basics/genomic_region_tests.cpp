// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include <boost/test/unit_test.hpp>

#include "basics/genomic_region.hpp"

namespace octopus { namespace test {

BOOST_AUTO_TEST_SUITE(basics)
BOOST_AUTO_TEST_SUITE(genomic_region)

BOOST_AUTO_TEST_CASE(comparing_order_of_regions_on_different_contigs_is_an_error)
{
    const GenomicRegion r1 {"1", 0, 1}, r2 {"2", 0, 1};
    BOOST_CHECK_NO_THROW(r1 == r2);
    BOOST_CHECK_NO_THROW(r1 != r2);
    BOOST_CHECK_THROW(r1 <= r2, RegionError);
    BOOST_CHECK_THROW(r1 < r2, RegionError);
    BOOST_CHECK_THROW(r1 > r2, RegionError);
    BOOST_CHECK_THROW(r1 >= r2, RegionError);
    BOOST_CHECK_THROW(is_before(r1, r2), RegionError);
    BOOST_CHECK_THROW(is_after(r1, r2), RegionError);
    BOOST_CHECK_THROW(begins_before(r1, r2), RegionError);
    BOOST_CHECK_THROW(ends_before(r1, r2), RegionError);
}

BOOST_AUTO_TEST_CASE(some_operations_are_well_defined_on_different_contigs)
{
    const GenomicRegion r1 {"1", 0, 1}, r2 {"2", 0, 1};
    BOOST_CHECK_NO_THROW(overlaps(r1, r2));
    BOOST_CHECK_NO_THROW(contains(r1, r2));
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
    
} // namespace test
} // namespace octopus
