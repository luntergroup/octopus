// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include <boost/test/unit_test.hpp>

#include <vector>

#include "basics/contig_region.hpp"
#include "concepts/mappable_range.hpp"

namespace octopus { namespace test {

BOOST_AUTO_TEST_SUITE(concepts)
BOOST_AUTO_TEST_SUITE(mappable_range)

BOOST_AUTO_TEST_CASE(mappable_ranges_can_be_empty)
{
    const std::vector<ContigRegion> v {};
    const ContigRegion r {};
    BOOST_REQUIRE_NO_THROW(make_overlap_range(v.cbegin(), v.cend(), r));
    const auto overlapped = make_overlap_range(v.cbegin(), v.cend(), r);
    BOOST_CHECK(overlapped.empty());
    BOOST_REQUIRE_NO_THROW(make_contained_range(v.cbegin(), v.cend(), r));
    const auto contained = make_contained_range(v.cbegin(), v.cend(), r);
    BOOST_CHECK(contained.empty());
}

BOOST_AUTO_TEST_CASE(single_element_ranges_are_ok)
{
    const std::vector<ContigRegion> v {ContigRegion {0, 1}};
    BOOST_REQUIRE_NO_THROW(make_overlap_range(v.cbegin(), v.cend(), v.front()));
    const auto overlapped = make_overlap_range(v.cbegin(), v.cend(), v.front());
    BOOST_CHECK(!overlapped.empty());
    BOOST_CHECK_EQUAL(overlapped.front(), v.front());
    BOOST_CHECK_EQUAL(overlapped.back(), v.back());
    BOOST_REQUIRE_NO_THROW(make_contained_range(v.cbegin(), v.cend(), v.front()));
    const auto contained = make_contained_range(v.cbegin(), v.cend(), v.front());
    BOOST_CHECK(!contained.empty());
    BOOST_CHECK_EQUAL(contained.front(), v.front());
    BOOST_CHECK_EQUAL(contained.back(), v.back());
}

BOOST_AUTO_TEST_CASE(overlap_ranges_filters_non_overlapped_elements)
{
    const std::vector<ContigRegion> v {
            ContigRegion {0, 5}, ContigRegion {1, 2}, ContigRegion {3, 4}
    };
    auto overlapped = make_overlap_range(v.cbegin(), v.cend(), v.back());
    BOOST_CHECK_EQUAL(overlapped.front(), v.front());
    BOOST_REQUIRE_NO_THROW(overlapped.drop_front());
    BOOST_CHECK_EQUAL(overlapped.front(), v.back());
    BOOST_REQUIRE_NO_THROW(overlapped.drop_front());
    BOOST_CHECK(overlapped.empty());
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

} // namespace test
} // namespace octopus
