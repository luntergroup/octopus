// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include <boost/test/unit_test.hpp>

#include <string>

#include "core/tools/vargen/utils/global_aligner.hpp"

namespace octopus { namespace test {
    
BOOST_AUTO_TEST_SUITE(core)
BOOST_AUTO_TEST_SUITE(global_aligner)

BOOST_AUTO_TEST_CASE(align_handles_empty_sequences)
{
    using coretools::align;
    
    constexpr coretools::Model defaultModel {};
    const std::string empty {}, nonempty {"A"};
    
    BOOST_REQUIRE_NO_THROW(align(empty, empty));
    auto alignment = align(empty, empty);
    BOOST_CHECK_EQUAL(alignment.cigar, CigarString {});
    BOOST_CHECK_EQUAL(alignment.score, 0);
    BOOST_REQUIRE_NO_THROW(align(empty, nonempty));
    
    alignment = align(empty, nonempty);
    BOOST_CHECK_EQUAL(alignment.cigar, parse_cigar("1I"));
    BOOST_CHECK_EQUAL(alignment.score, defaultModel.gap_open);
    BOOST_REQUIRE_NO_THROW(align(nonempty, empty));
    
    alignment = align(nonempty, empty);
    BOOST_CHECK_EQUAL(alignment.cigar, parse_cigar("1D"));
    BOOST_CHECK_EQUAL(alignment.score, defaultModel.gap_open);
}

BOOST_AUTO_TEST_CASE(align_returns_the_optimal_global_alignmenet)
{
    using coretools::align;
    
    constexpr coretools::Model defaultModel {};
    const std::string sequence1 {"AAA"}, sequence2 {"ACA"};
    
    auto alignment = align(sequence1, sequence2);
    BOOST_CHECK_EQUAL(alignment.cigar, parse_cigar("1=1X1="));
    BOOST_CHECK_EQUAL(alignment.score, 2 * defaultModel.match + defaultModel.mismatch);
    
    const std::string sequence3 {"AAAA"}, sequence4 {"AA"};
    alignment = align(sequence1, sequence3);
    BOOST_CHECK_EQUAL(alignment.cigar, parse_cigar("1I3="));
    BOOST_CHECK_EQUAL(alignment.score, 3 * defaultModel.match + defaultModel.gap_open);
    
    alignment = align(sequence1, sequence4);
    BOOST_CHECK_EQUAL(alignment.cigar, parse_cigar("1D2="));
    BOOST_CHECK_EQUAL(alignment.score, 2 * defaultModel.match + defaultModel.gap_open);
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
    
} // namespace test
} // namespace octopus
