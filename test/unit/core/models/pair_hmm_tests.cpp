// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include <boost/test/unit_test.hpp>

#include <string>
#include <vector>
#include <cstddef>
#include <algorithm>
#include <utility>
#include <iostream>

#include "core/models/pairhmm/simd_pair_hmm.hpp"

namespace octopus { namespace test {

using namespace octopus::hmm::simd;

BOOST_AUTO_TEST_SUITE(core)
BOOST_AUTO_TEST_SUITE(model)

struct TestCase
{
    std::string target, query;
    std::vector<std::int8_t> base_qualities;
    std::vector<std::int8_t> gap_open;
    int gap_extend, nuc_prior;
    unsigned band_size = 8;
};

struct Alignment
{
    int score;
    int begin;
    std::string target, query;
};

std::ostream& operator<<(std::ostream& os, const std::vector<std::int8_t>& quals)
{
    std::copy(std::cbegin(quals), std::cend(quals), std::ostream_iterator<std::int8_t> {os, " "});
    return os;
}

std::ostream& operator<<(std::ostream& os, const TestCase& test)
{
    os << "target: " << test.target << '\n';
    os << "query: " << test.query << '\n';
    os << "base_qualities: " << test.base_qualities << '\n';
    os << "gap_open: " << test.gap_open << '\n';
    os << "gap_extend: " << test.gap_extend << '\n';
    os << "nuc_prior: " << test.nuc_prior << '\n';
    return os;
}

std::ostream& operator<<(std::ostream& os, const Alignment& alignment)
{
    os << "score: " << alignment.score << '\n';
    os << "begin: " << alignment.begin << '\n';
    os << "target: " << alignment.target << '\n';
    os << "query: " << alignment.query << '\n';
    return os;
}

template <typename SIMD>
auto
align_score_helper(TestCase test)
{
    return align<SIMD>(test.target.data(), test.query.data(), test.base_qualities.data(),
                       static_cast<int>(test.target.size()), static_cast<int>(test.query.size()),
                       test.gap_open.data(), test.gap_extend, test.nuc_prior);
}

template <typename SIMD>
Alignment
align_helper(TestCase test)
{
    Alignment result {};
    std::vector<char> align1(2 * test.target.size() + 1, '\0'), align2(2 * test.target.size() + 1, '\0');
    result.score = align<SIMD>(test.target.data(), test.query.data(), test.base_qualities.data(),
                               static_cast<int>(test.target.size()), static_cast<int>(test.query.size()),
                               test.gap_open.data(), test.gap_extend, test.nuc_prior,
                               result.begin, align1.data(), align2.data());
    result.target.assign(align1.cbegin(), std::find(align1.cbegin(), align1.cend(), '\0'));
    result.query.assign(align2.cbegin(), std::find(align2.cbegin(), align2.cend(), '\0'));
    return result;
}

#define CHECK_ALIGNMENT(received, expected) \
    BOOST_CHECK_EQUAL(received.score, expected.score); \
    BOOST_CHECK_EQUAL(received.begin, expected.begin); \
    BOOST_CHECK_EQUAL(received.target, expected.target); \
    BOOST_CHECK_EQUAL(received.query, expected.query);

BOOST_AUTO_TEST_CASE(sse2_check_alignments)
{
    TestCase test;
    Alignment expected_alignment;
    
    // test 1
    test = {
        "ACGTACGTACGTACGAAAA",
        "AAAA",
        {40,40,40,40},
        {10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10},
        1,
        4
    };
    expected_alignment = {
        0,
        15,
        "AAAA",
        "AAAA"
    };
    BOOST_CHECK_EQUAL(align_score_helper<SSE2>(test), expected_alignment.score);
    CHECK_ALIGNMENT(align_helper<SSE2>(test), expected_alignment);
    
    // test 2
    test = {
        "ACGTACGTACGTACGAATA",
        "AAAA",
        {40,40,40,40},
        {90,90,90,90,90,90,90,90,90,90,90,90,90,90,90,90,90,90,90},
        1,
        4
    };
    expected_alignment = {
        40,
        15,
        "AATA",
        "AAAA"
    };
    BOOST_CHECK_EQUAL(align_score_helper<SSE2>(test), expected_alignment.score);
    CHECK_ALIGNMENT(align_helper<SSE2>(test), expected_alignment);
    
    // test 3
    test = {
        "ACGTACGTACGTACGAAGC",
        "CGGC",
        {40,40,40,40},
        {90,90,90,90,90,90,90,90,90,90,90,90,90,90,90,70,90,90,90},
        1,
        4
    };
    expected_alignment = {
        71,
        13,
        "CGAAGC",
        "CG--GC"
    };
    BOOST_CHECK_EQUAL(align_score_helper<SSE2>(test), expected_alignment.score);
    CHECK_ALIGNMENT(align_helper<SSE2>(test), expected_alignment);
    
    // test 4
    test = {
        "CGAAGCACGTACGTACGTA",
        "CGGC",
        {40,40,40,40},
        {90,90,70,90,90,90,90,90,90,90,90,90,90,90,90,90,90,90,90},
        1,
        4
    };
    expected_alignment = {
        71,
        0,
        "CGAAGC",
        "CG--GC"
    };
    BOOST_CHECK_EQUAL(align_score_helper<SSE2>(test), expected_alignment.score);
    CHECK_ALIGNMENT(align_helper<SSE2>(test), expected_alignment);
}

#ifdef __AVX2__
BOOST_AUTO_TEST_CASE(avx2_check_alignments)
{
    TestCase test;
    Alignment expected_alignment;
    
    // test 1
    test = {
        "ACGTACGTACGTACGTACGTACGTACGTACGAAAA",
        "AAAA",
        {40,40,40,40},
        {10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10,10},
        1,
        4
    };
    expected_alignment = {
        0,
        15,
        "AAAA",
        "AAAA"
    };
    BOOST_CHECK_EQUAL(align_score_helper<AVX2>(test), expected_alignment.score);
    CHECK_ALIGNMENT(align_helper<AVX2>(test), expected_alignment);
    
    // test 2
    test = {
        "ACGTACGTACGTACGTACGTACGTACGTACGAATA",
        "AAAA",
        {40,40,40,40},
        {90,90,90,90,90,90,90,90,90,90,90,90,90,90,90,90,90,90,90},
        1,
        4
    };
    expected_alignment = {
        40,
        15,
        "AATA",
        "AAAA"
    };
    BOOST_CHECK_EQUAL(align_score_helper<AVX2>(test), expected_alignment.score);
    CHECK_ALIGNMENT(align_helper<AVX2>(test), expected_alignment);
    
    // test 3
    test = {
        "ACGTACGTACGTACGAAGCACGTACGTACGTACGT",
        "CGGC",
        {40,40,40,40},
        {90,90,90,90,90,90,90,90,90,90,90,90,90,90,90,70,90,90,90},
        1,
        4
    };
    expected_alignment = {
        71,
        13,
        "CGAAGC",
        "CG--GC"
    };
    BOOST_CHECK_EQUAL(align_score_helper<AVX2>(test), expected_alignment.score);
    CHECK_ALIGNMENT(align_helper<AVX2>(test), expected_alignment);
    
    // test 4
    test = {
        "CGAAGCACGTACGTACGTAACGTACGTACGTACGT",
        "CGGC",
        {40,40,40,40},
        {90,90,70,90,90,90,90,90,90,90,90,90,90,90,90,90,90,90,90},
        1,
        4
    };
    expected_alignment = {
        71,
        0,
        "CGAAGC",
        "CG--GC"
    };
    BOOST_CHECK_EQUAL(align_score_helper<AVX2>(test), expected_alignment.score);
    CHECK_ALIGNMENT(align_helper<AVX2>(test), expected_alignment);
}
#endif /* __AVX2__ */

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

} // namespace test
} // namespace octopus
