// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include <boost/test/unit_test.hpp>

#include <string>

#include "basics/cigar_string.hpp"

namespace octopus { namespace test {

BOOST_AUTO_TEST_SUITE(basics)
BOOST_AUTO_TEST_SUITE(cigar_string)

BOOST_AUTO_TEST_CASE(cigars_with_the_same_ordered_ops_are_equal)
{
    using CO   = CigarOperation;
    using Flag = CO::Flag;
    
    CigarString cigar1 {}, cigar2 {CO {10, Flag::alignmentMatch}}, cigar3 {CO {10, Flag::sequenceMatch}};
    
    BOOST_CHECK_EQUAL(cigar1, cigar1);
    BOOST_CHECK_EQUAL(cigar2, cigar2);
    BOOST_CHECK_EQUAL(cigar3, cigar3);
    
    BOOST_CHECK_NE(cigar1, cigar2);
    BOOST_CHECK_NE(cigar1, cigar3);
    BOOST_CHECK_NE(cigar2, cigar3);
}

BOOST_AUTO_TEST_CASE(parse_cigar_works)
{
    using CO   = CigarOperation;
    using Flag = CO::Flag;
    
    CigarString cigar1 {}, cigar2 {CO {10, Flag::alignmentMatch}}, cigar3 {CO {5, Flag::insertion}};
    BOOST_CHECK_NO_THROW(parse_cigar(""));
    
    auto pcigar1 = parse_cigar("");
    BOOST_CHECK_EQUAL(pcigar1, cigar1);
    
    BOOST_CHECK_NO_THROW(parse_cigar("10M"));
    auto pcigar2 = parse_cigar("10M");
    BOOST_CHECK_EQUAL(pcigar2, cigar2);
    
    auto pcigar3 = parse_cigar("5I");
    BOOST_CHECK_EQUAL(pcigar3, cigar3);
}

BOOST_AUTO_TEST_CASE(a_cigar_is_valid_if_all_ops_are_valid)
{
    auto cigar1 = parse_cigar("");
    BOOST_CHECK(!is_valid(cigar1));
    
    auto cigar2 = parse_cigar("10M");
    auto cigar3 = parse_cigar("5I10M");
    auto cigar4 = parse_cigar("10M10M");
    BOOST_CHECK(is_valid(cigar2));
    BOOST_CHECK(is_valid(cigar3));
    BOOST_CHECK(is_valid(cigar4));
    
    auto cigar5 = parse_cigar("5S1D19M9I2I4D28X1=1D6S10H");
    BOOST_CHECK(is_valid(cigar5));
    
    auto cigar6 = parse_cigar("1T");
    BOOST_CHECK(!is_valid(cigar6));
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
    
} // namespace test
} // namespace octopus
