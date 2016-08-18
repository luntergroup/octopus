// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include <boost/test/unit_test.hpp>

#include <vector>
#include <iterator>
#include <algorithm>
#include <future>

#include <io/region/region_parser.hpp>
#include <exceptions/user_error.hpp>

#include "mock/mock_reference.hpp"

namespace octopus { namespace test {

BOOST_AUTO_TEST_SUITE(io)
BOOST_AUTO_TEST_SUITE(region_parser)

BOOST_AUTO_TEST_CASE(parse_region_throws_when_given_badly_formated_regions)
{
    using ::octopus::io::parse_region;
    
    const auto reference = mock::make_reference();
    
    BOOST_CHECK_THROW(parse_region("", reference), UserError);
    BOOST_CHECK_THROW(parse_region("-", reference), UserError);
    BOOST_CHECK_THROW(parse_region("5:100-99", reference), UserError);
    BOOST_CHECK_THROW(parse_region("not_in_reference", reference), UserError);
    BOOST_CHECK_THROW(parse_region("0", reference), UserError);
    BOOST_CHECK_THROW(parse_region("-1", reference), UserError);
    BOOST_CHECK_THROW(parse_region("--1", reference), UserError);
    BOOST_CHECK_THROW(parse_region("1:", reference), UserError);
    BOOST_CHECK_THROW(parse_region("1:-", reference), UserError);
    BOOST_CHECK_THROW(parse_region("1:-0-10", reference), UserError);
    BOOST_CHECK_THROW(parse_region("1:-343-1000", reference), UserError);
    BOOST_CHECK_THROW(parse_region("1:0--10", reference), UserError);
    BOOST_CHECK_THROW(parse_region("1:0--010", reference), UserError);
    BOOST_CHECK_THROW(parse_region("1:o-1", reference), UserError);
    BOOST_CHECK_THROW(parse_region("1:0-1o", reference), UserError);
    BOOST_CHECK_THROW(parse_region("1:0-1o0", reference), UserError);
    BOOST_CHECK_THROW(parse_region("2::0-323", reference), UserError);
    BOOST_CHECK_THROW(parse_region("3:0:-1", reference), UserError);
    BOOST_CHECK_THROW(parse_region("3:0:-948", reference), UserError);
    BOOST_CHECK_THROW(parse_region("3:0-:10", reference), UserError);
    BOOST_CHECK_THROW(parse_region("3:0--10", reference), UserError);
    BOOST_CHECK_THROW(parse_region("3:-0-10", reference), UserError);
    BOOST_CHECK_THROW(parse_region("3:-1-10", reference), UserError);
    BOOST_CHECK_THROW(parse_region("3:o:-1o0", reference), UserError);
    BOOST_CHECK_THROW(parse_region("3:0:-1o0", reference), UserError);
    BOOST_CHECK_THROW(parse_region("3:0:-1o0", reference), UserError);
}

BOOST_AUTO_TEST_CASE(parse_region_works_with_correctly_formatted_region_input)
{
    using ::octopus::io::parse_region;
    
    const auto reference = mock::make_reference();
    
    GenomicRegion region;
    
    BOOST_REQUIRE_NO_THROW(region = parse_region("1", reference));
    
    BOOST_CHECK(region.contig_name() == "1");
    BOOST_CHECK(region.begin() == 0);
    BOOST_CHECK(region.end() == reference.contig_size("1"));
    
    BOOST_REQUIRE_NO_THROW(region = parse_region("1:100-200", reference));
    
    BOOST_CHECK(region.contig_name() == "2");
    BOOST_CHECK(region.begin() == 100);
    BOOST_CHECK(region.end() == 200);
    
    BOOST_REQUIRE_NO_THROW(region = parse_region("4:1,121-2,491", reference));
    
    BOOST_CHECK(region.contig_name() == "4");
    BOOST_CHECK(region.begin() == 1121);
    BOOST_CHECK(region.end() == 2491);
    
    BOOST_REQUIRE_NO_THROW(region = parse_region("3:99-", reference));
    
    BOOST_CHECK(region.contig_name() == "3");
    BOOST_CHECK(region.begin() == 99);
    BOOST_CHECK(region.end() == reference.contig_size("3"));
    
    BOOST_REQUIRE_NO_THROW(region = parse_region("5:3", reference));
    
    BOOST_CHECK(region.contig_name() == "5");
    BOOST_CHECK(region.begin() == 3);
    BOOST_CHECK(region.end() == 4);
    
    BOOST_REQUIRE_NO_THROW(region = parse_region("6:00-0100", reference));
    
    BOOST_CHECK(region.contig_name() == "6");
    BOOST_CHECK(region.begin() == 0);
    BOOST_CHECK(region.end() == 100);
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

} // namespace test
} // namespace octopus
