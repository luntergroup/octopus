// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include <boost/test/unit_test.hpp>

#include <utility>

#include "basics/contig_region.hpp"
#include "core/types/allele.hpp"

namespace octopus { namespace test {

template <typename T>
auto make_allele(ContigRegion::Position begin, ContigRegion::Position end, T&& sequence)
{
    return ContigAllele {ContigRegion {begin, end}, std::forward<T>(sequence)};
}

BOOST_AUTO_TEST_SUITE(core)
BOOST_AUTO_TEST_SUITE(allele)

BOOST_AUTO_TEST_CASE(alleles_are_equal_if_their_region_and_sequence_match)
{
    auto a1 = make_allele(0, 1, "A"), a2 = make_allele(0, 1, "C"), a3 = make_allele(0, 0, ""),
         a4 = make_allele(0, 0, "A"), a5 = make_allele(0, 1, "AA"), a6 = make_allele(0, 1, "AC"),
         a7 = make_allele(0, 1, "");
    
    BOOST_CHECK_EQUAL(a1, a1);
    BOOST_CHECK_EQUAL(a2, a2);
    BOOST_CHECK_EQUAL(a3, a3);
    BOOST_CHECK_EQUAL(a4, a4);
    BOOST_CHECK_EQUAL(a5, a5);
    BOOST_CHECK_EQUAL(a6, a6);
    BOOST_CHECK_EQUAL(a7, a7);
    
    BOOST_CHECK_NE(a1, a2);
    BOOST_CHECK_NE(a1, a3);
    BOOST_CHECK_NE(a2, a3);
    BOOST_CHECK_NE(a1, a4);
    BOOST_CHECK_NE(a1, a5);
    BOOST_CHECK_NE(a1, a6);
    BOOST_CHECK_NE(a1, a7);
    BOOST_CHECK_NE(a2, a4);
    BOOST_CHECK_NE(a2, a5);
    BOOST_CHECK_NE(a2, a6);
    BOOST_CHECK_NE(a2, a7);
    BOOST_CHECK_NE(a3, a2);
    BOOST_CHECK_NE(a3, a4);
    BOOST_CHECK_NE(a3, a5);
    BOOST_CHECK_NE(a3, a6);
    BOOST_CHECK_NE(a3, a7);
    BOOST_CHECK_NE(a4, a5);
    BOOST_CHECK_NE(a4, a6);
    BOOST_CHECK_NE(a4, a7);
    BOOST_CHECK_NE(a5, a6);
    BOOST_CHECK_NE(a6, a7);
}

BOOST_AUTO_TEST_CASE(operators_less_and_equal_are_consistent)
{
    
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

} // namespace test
} // namespace octopus
