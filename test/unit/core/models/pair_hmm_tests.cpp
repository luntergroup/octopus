// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include <boost/test/unit_test.hpp>

#include <string>
#include <vector>
#include <cstddef>
#include <algorithm>
#include <utility>
#include <iostream>

#include "core/models/pairhmm/simd_pair_hmm_factory.hpp"

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

template <typename HMM>
auto
align_score_helper(TestCase test, HMM hmm)
{
    return hmm.align(test.target.data(), test.query.data(), test.base_qualities.data(),
                     static_cast<int>(test.target.size()), static_cast<int>(test.query.size()),
                     test.gap_open.data(), test.gap_extend, test.nuc_prior);
}

template <typename HMM>
Alignment
align_helper(TestCase test, HMM hmm)
{
    Alignment result {};
    std::vector<char> align1(2 * test.target.size() + 1, '\0'), align2(2 * test.target.size() + 1, '\0');
    result.score = hmm.align(test.target.data(), test.query.data(), test.base_qualities.data(),
                             static_cast<int>(test.target.size()), static_cast<int>(test.query.size()),
                             test.gap_open.data(), test.gap_extend, test.nuc_prior,
                             result.begin, align1.data(), align2.data());
    result.target.assign(align1.cbegin(), std::find(align1.cbegin(), align1.cend(), '\0'));
    result.query.assign(align2.cbegin(), std::find(align2.cbegin(), align2.cend(), '\0'));
    return result;
}

#define CHECK_TEST(test, hmm) \
        BOOST_CHECK_EQUAL(test.target.size(), test.query.size() + 2 * hmm.band_size() - 1); \
        BOOST_CHECK_EQUAL(test.query.size(), test.base_qualities.size()); \
        BOOST_CHECK_EQUAL(test.target.size(), test.gap_open.size());

#define CHECK_ALIGNER(test, hmm, expected_alignment) \
        BOOST_CHECK_EQUAL(align_score_helper(test, hmm), expected_alignment.score); \
        BOOST_CHECK_EQUAL(align_helper(test, hmm).score, expected_alignment.score); \
        BOOST_CHECK_EQUAL(align_helper(test, hmm).begin, expected_alignment.begin); \
        BOOST_CHECK_EQUAL(align_helper(test, hmm).target, expected_alignment.target); \
        BOOST_CHECK_EQUAL(align_helper(test, hmm).query, expected_alignment.query);

#define CHECK_SPEED(test, hmm, expected_alignment, iterations)     \
        CHECK_TEST(test, hmm) \
        CHECK_ALIGNER(test, hmm, expected_alignment) \
        for(int i=0; i<iterations; i++) { BOOST_CHECK_EQUAL(align_score_helper(test, hmm), expected_alignment.score); }

TestCase band8_speed_test = {
    "AACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAACCCCACGTATATATATATATATGGGGACGTAACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAA",
    "AACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAACCCCACGTGGGACGTAACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAA",
    {40,40,40,40, 40,40,40,40, 40,40,40,40, 40,40,40,40,
     40,40,40,40, 40,40,40,40, 40,40,40,40, 40,40,40,40,
     40,40,40,40, 40,40,40,40, 40,40,40,40, 40,40,40,40,
     40,40,40,40, 40,40,40,40, 40,40,40,40, 40,40,40,40,
     40,40,40,40, 40,40,40,40, 40,40,40,40, 40,40,40,40,
     40,40,40,40, 40,40,40,40, 40,40,40,40, 40,40,40,40,
     40,40,40,40, 40,40,40,40, 40,40,40,40, 40,40,40,40,
     40,40,40,40, 40,40,40,40, 40,40,40,40, 40,40,40,40,
     40,40,40,40, 40,40,40,40, 40,40,40,40, 40,40,40,40,
     40,40,40,40, 40,40,40,40, 40,40,40,40, 40,40,40,40,
     40,40,40,40, 40,40,40,40, 40,40,40,40, 40,40,40,40,
     40,40,40,40, 40,40,40,40, 40,40,40,40, 40,40,40,40,
     40,40,40,40, 40,40,40,40, 40,40,40,40, 40,40,40,40,
     40,40,40,40, 40,40,40,40, 40,40,40,40, 40,40,40,40,
     40,40,40,40, 40,40,40,40, 40,40,40,40, 40,40,40,40,
     40,40,40,40, 40,40,40,40, 40,40,40,40, 40,40,40,40,
     40,40,40,40,40,40,40,40,40,40,40,40,40,40,40},
    {90,90,90,90, 90,90,90,90, 90,90,90,90, 90,90,90,90,
     90,90,90,90, 90,90,90,90, 90,90,90,90, 90,90,90,90,
     90,90,90,90, 90,90,90,90, 90,90,90,90, 90,90,90,90,
     90,90,90,90, 90,90,90,90, 90,90,90,90, 90,90,90,90,
     90,90,90,90, 90,90,90,90, 90,90,90,90, 90,90,90,90,
     90,90,90,90, 90,90,90,90, 90,90,90,90, 90,90,90,90,
     90,90,90,90, 90,90,90,90, 90,90,90,90, 90,90,90,90,
     90,90,90,90, 90,90,90,90, 90,90,90,90, 90,90,90,90,
     90,90,90,90, 90,90,90,90, 70,90,90,90, 90,90,90,90,
     90,90,90,90, 90,90,90,90, 90,90,90,90, 90,90,90,90,
     90,90,90,90, 90,90,90,90, 90,90,90,90, 90,90,90,90,
     90,90,90,90, 90,90,90,90, 90,90,90,90, 90,90,90,90,
     90,90,90,90, 90,90,90,90, 90,90,90,90, 90,90,90,90,
     90,90,90,90, 90,90,90,90, 90,90,90,90, 90,90,90,90,
     90,90,90,90, 90,90,90,90, 90,90,90,90, 90,90,90,90,
     90,90,90,90, 90,90,90,90, 90,90,90,90, 90,90,90,90,
     90,90,90,90, 90,90,90,90, 90,90,90,90, 90,90,90,90,
     90,90,90,90, 90,90,90,90, 90,90,90,90, 90,90},
    1,
    4
};
Alignment band8_speed_expected_alignment = {
    84,
    0,
    "AACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAACCCCACGTATATATATATATATGGGGACGTAACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAA",
    "AACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAACCCCACGT---------------GGGACGTAACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAA"
};

TestCase band16_speed_test = {
    "AACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAACCCCACGTATATATATATATATGGGGACGTAACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAA",
    "AACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAACCCCACGTGGGACGTAACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAA",
    {40,40,40,40, 40,40,40,40, 40,40,40,40, 40,40,40,40,
     40,40,40,40, 40,40,40,40, 40,40,40,40, 40,40,40,40,
     40,40,40,40, 40,40,40,40, 40,40,40,40, 40,40,40,40,
     40,40,40,40, 40,40,40,40, 40,40,40,40, 40,40,40,40,
     40,40,40,40, 40,40,40,40, 40,40,40,40, 40,40,40,40,
     40,40,40,40, 40,40,40,40, 40,40,40,40, 40,40,40,40,
     40,40,40,40, 40,40,40,40, 40,40,40,40, 40,40,40,40,
     40,40,40,40, 40,40,40,40, 40,40,40,40, 40,40,40,40,
     40,40,40,40, 40,40,40,40, 40,40,40,40, 40,40,40,40,
     40,40,40,40, 40,40,40,40, 40,40,40,40, 40,40,40,40,
     40,40,40,40, 40,40,40,40, 40,40,40,40, 40,40,40,40,
     40,40,40,40, 40,40,40,40, 40,40,40,40, 40,40,40,40,
     40,40,40,40, 40,40,40,40, 40,40,40,40, 40,40,40,40,
     40,40,40,40, 40,40,40,40, 40,40,40,40, 40,40,40,40,
     40,40,40,40, 40,40,40,40, 40,40,40,40, 40,40,40,40,
     40,40,40,40, 40,40,40,40, 40,40,40,40, 40,40,40,40,
     40,40,40,40,40,40,40,40,40,40,40,40,40,40,40},
    {90,90,90,90, 90,90,90,90, 90,90,90,90, 90,90,90,90,
     90,90,90,90, 90,90,90,90, 90,90,90,90, 90,90,90,90,
     90,90,90,90, 90,90,90,90, 90,90,90,90, 90,90,90,90,
     90,90,90,90, 90,90,90,90, 90,90,90,90, 90,90,90,90,
     90,90,90,90, 90,90,90,90, 90,90,90,90, 90,90,90,90,
     90,90,90,90, 90,90,90,90, 90,90,90,90, 90,90,90,90,
     90,90,90,90, 90,90,90,90, 90,90,90,90, 90,90,90,90,
     90,90,90,90, 90,90,90,90, 90,90,90,90, 90,90,90,90,
     90,90,90,90, 90,90,90,90, 70,90,90,90, 90,90,90,90,
     90,90,90,90, 90,90,90,90, 90,90,90,90, 90,90,90,90,
     90,90,90,90, 90,90,90,90, 90,90,90,90, 90,90,90,90,
     90,90,90,90, 90,90,90,90, 90,90,90,90, 90,90,90,90,
     90,90,90,90, 90,90,90,90, 90,90,90,90, 90,90,90,90,
     90,90,90,90, 90,90,90,90, 90,90,90,90, 90,90,90,90,
     90,90,90,90, 90,90,90,90, 90,90,90,90, 90,90,90,90,
     90,90,90,90, 90,90,90,90, 90,90,90,90, 90,90,90,90,
     90,90,90,90, 90,90,90,90, 90,90,90,90, 90,90,90,90,
     90,90,90,90, 90,90,90,90, 90,90,90,90, 90,90,90,90,
     90,90,90,90, 90,90,90,90, 90,90,90,90, 90,90},
    1,
    4
};
Alignment band16_speed_expected_alignment = {
    84,
    0,
    "AACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAACCCCACGTATATATATATATATGGGGACGTAACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAA",
    "AACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAACCCCACGT---------------GGGACGTAACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAAAACCGGTTTTGGCCAA"
};



BOOST_AUTO_TEST_CASE(sse2_band_size_8_check_alignments)
{
    SSE2PairHMM<8, short> sse_short_hmm;
    SSE2PairHMM<8, int> sse_int_hmm;
    TestCase test;
    Alignment expected_alignment;
    
    // test 1
    test = {
        "ACGTACGTACGTACGAAAA",
        "AAAA",
        {40,40,40,40},
        {10,10,10,10, 10,10,10,10, 10,10,10,10, 10,10,10,10,
         10,10,10},
        1,
        4
    };
    expected_alignment = {
        0,
        15,
        "AAAA",
        "AAAA"
    };
    CHECK_TEST(test, sse_short_hmm)
    CHECK_ALIGNER(test, sse_short_hmm, expected_alignment)
    CHECK_TEST(test, sse_int_hmm)
    CHECK_ALIGNER(test, sse_int_hmm, expected_alignment)
    
    // test 2
    test = {
        "ACGTACGTACGTACGAATA",
        "AAAA",
        {40,40,40,40},
        {90,90,90,90, 90,90,90,90, 90,90,90,90, 90,90,90,90,
         90,90,90},
        1,
        4
    };
    expected_alignment = {
        40,
        15,
        "AATA",
        "AAAA"
    };
    CHECK_TEST(test, sse_short_hmm)
    CHECK_ALIGNER(test, sse_short_hmm, expected_alignment)
    CHECK_TEST(test, sse_int_hmm)
    CHECK_ALIGNER(test, sse_int_hmm, expected_alignment)
    
    // test 3
    test = {
        "ACGTACGAAGCTACGTACG",
        "CGGC",
        {40,40,40,40},
        {90,90,90,90, 90,90,90,70, 90,90,90,90, 90,90,90,90,
         90,90,90},
        1,
        4
    };
    expected_alignment = {
        71,
        5,
        "CGAAGC",
        "CG--GC"
    };
    CHECK_TEST(test, sse_short_hmm)
    CHECK_ALIGNER(test, sse_short_hmm, expected_alignment)
    CHECK_TEST(test, sse_int_hmm)
    CHECK_ALIGNER(test, sse_int_hmm, expected_alignment)
    
    // test 4
    test = {
        "CGAAGCACGTACGTACGTA",
        "CGGC",
        {40,40,40,40},
        {90,90,70,90, 90,90,90,90, 90,90,90,90, 90,90,90,90, 
         90,90,90},
        1,
        4
    };
    expected_alignment = {
        71,
        0,
        "CGAAGC",
        "CG--GC"
    };
    CHECK_TEST(test, sse_short_hmm)
    CHECK_ALIGNER(test, sse_short_hmm, expected_alignment)
    CHECK_TEST(test, sse_int_hmm)
    CHECK_ALIGNER(test, sse_int_hmm, expected_alignment)

    // test 5
    test = {
        "CCCCACGTATATATATATATATGGGGACGT",
        "CCCCACGTGGGACGT",
        {40,40,40,40,40,40,40,40,40,40,40,40,40,40,40},
        {90,90,90,90, 90,90,90,90, 70,90,90,90, 90,90,90,90,
         90,90,90,90, 90,90,90,90, 90,90,90,90, 90,90},
        1,
        4
    };
    expected_alignment = {
        84,
        0,
        "CCCCACGTATATATATATATATGGGGACGT",
        "CCCCACGT---------------GGGACGT"
    };
    CHECK_TEST(test, sse_short_hmm)
    CHECK_ALIGNER(test, sse_short_hmm, expected_alignment)
    CHECK_TEST(test, sse_int_hmm)
    CHECK_ALIGNER(test, sse_int_hmm, expected_alignment)
}


BOOST_AUTO_TEST_CASE(sse2_band_size_16_check_alignments)
{
    SSE2PairHMM<16, short> sse_short_hmm;
    SSE2PairHMM<16, int> sse_int_hmm;
    TestCase test;
    Alignment expected_alignment;
    
    // test 1
    test = {
        "ACGTACGTACGTACGTACGTACGTACGTACGAAAA",
        "AAAA",
        {40,40,40,40},
        {10,10,10,10,10,10,10,10,
         10,10,10,10,10,10,10,10,
         10,10,10,10,10,10,10,10,
         10,10,10,10,10,10,10,10,
         10,10,10},
        1,
        4
    };
    expected_alignment = {
        0,
        31,
        "AAAA",
        "AAAA"
    };
    CHECK_TEST(test, sse_short_hmm)
    CHECK_ALIGNER(test, sse_short_hmm, expected_alignment)
    CHECK_TEST(test, sse_int_hmm)
    CHECK_ALIGNER(test, sse_int_hmm, expected_alignment)
    
    // test 2
    test = {
        "ACGTACGTACGTACGTACGTACGTACGTACGAATA",
        "AAAA",
        {40,40,40,40},
        {90,90,90,90,90,90,90,90,
         90,90,90,90,90,90,90,90,
         90,90,90,90,90,90,90,90,
         90,90,90,90,90,90,90,90,
         90,90,90},
        1,
        4
    };
    expected_alignment = {
        40,
        31,
        "AATA",
        "AAAA"
    };
    CHECK_TEST(test, sse_short_hmm)
    CHECK_ALIGNER(test, sse_short_hmm, expected_alignment)
    CHECK_TEST(test, sse_int_hmm)
    CHECK_ALIGNER(test, sse_int_hmm, expected_alignment)
    
    // test 3
    test = {
        "ACGTACGTACGTACGAAGCACGTACGTACGTACGT",
        "CGGC",
        {40,40,40,40},
        {90,90,90,90,90,90,90,90,
         90,90,90,90,90,90,90,70,
         90,90,90,90,90,90,90,90,
         90,90,90,90,90,90,90,90,
         90,90,90},
        1,
        4
    };
    expected_alignment = {
        71,
        13,
        "CGAAGC",
        "CG--GC"
    };
    CHECK_TEST(test, sse_short_hmm)
    CHECK_ALIGNER(test, sse_short_hmm, expected_alignment)
    CHECK_TEST(test, sse_int_hmm)
    CHECK_ALIGNER(test, sse_int_hmm, expected_alignment)
    
    // test 4
    test = {
        "CGAAGCACGTACGTACGTAACGTACGTACGTACGT",
        "CGGC",
        {40,40,40,40},
        {90,90,70,90,90,90,90,90,
         90,90,90,90,90,90,90,90,
         90,90,90,90,90,90,90,90,
         90,90,90,90,90,90,90,90,
         90,90,90},
        1,
        4
    };
    expected_alignment = {
        71,
        0,
        "CGAAGC",
        "CG--GC"
    };
    CHECK_TEST(test, sse_short_hmm)
    CHECK_ALIGNER(test, sse_short_hmm, expected_alignment)
    CHECK_TEST(test, sse_int_hmm)
    CHECK_ALIGNER(test, sse_int_hmm, expected_alignment)
    
    // test 5
    test = {
        "CCCCACGTCCCCACGTATATATATATATATGGGGACGTGGGGACGT",
        "CCCCACGTGGGACGT",
        {40,40,40,40,40,40,40,40,40,40,40,40,40,40,40},
        {90,90,90,90,90,90,90,90,
         90,90,90,90,90,90,90,90,
         70,90,90,90,90,90,90,90,
         90,90,90,90,90,90,90,90,
         90,90,90,90,90,90,90,90,
         90,90,90,90,90,90},
        1,
        4
    };
    expected_alignment = {
        84,
        8,
        "CCCCACGTATATATATATATATGGGGACGT",
        "CCCCACGT---------------GGGACGT"
    };
    CHECK_TEST(test, sse_short_hmm)
    CHECK_ALIGNER(test, sse_short_hmm, expected_alignment)
    CHECK_TEST(test, sse_int_hmm)
    CHECK_ALIGNER(test, sse_int_hmm, expected_alignment)
}

BOOST_AUTO_TEST_CASE(sse2_band_size_32_check_alignments)
{
    SSE2PairHMM<32, short> sse_short_hmm;
    SSE2PairHMM<32, int> sse_int_hmm;
    TestCase test;
    Alignment expected_alignment;
    
    // test 1
    test = {
        "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGAAAA",
        "AAAA",
        {40,40,40,40},
        {10,10,10,10,10,10,10,10,
         10,10,10,10,10,10,10,10,
         10,10,10,10,10,10,10,10,
         10,10,10,10,10,10,10,10,
         10,10,10,10,10,10,10,10,
         10,10,10,10,10,10,10,10,
         10,10,10,10,10,10,10,10,
         10,10,10,10,10,10,10,10,
         10,10,10},
        1,
        4
    };
    expected_alignment = {
        0,
        63,
        "AAAA",
        "AAAA"
    };
    CHECK_TEST(test, sse_short_hmm)
    CHECK_ALIGNER(test, sse_short_hmm, expected_alignment)
    CHECK_TEST(test, sse_int_hmm)
    CHECK_ALIGNER(test, sse_int_hmm, expected_alignment)

    // test 2
    test = {
        "ACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGAATA",
        "AAAA",
        {40,40,40,40},
        {90,90,90,90,90,90,90,90,
         90,90,90,90,90,90,90,90,
         90,90,90,90,90,90,90,90,
         90,90,90,90,90,90,90,90,
         90,90,90,90,90,90,90,90,
         90,90,90,90,90,90,90,90,
         90,90,90,90,90,90,90,90,
         90,90,90,90,90,90,90,90,
         90,90,90},
        1,
        4
    };
    expected_alignment = {
        40,
        63,
        "AATA",
        "AAAA"
    };
    CHECK_TEST(test, sse_short_hmm)
    CHECK_ALIGNER(test, sse_short_hmm, expected_alignment)
    CHECK_TEST(test, sse_int_hmm)
    CHECK_ALIGNER(test, sse_int_hmm, expected_alignment)
    
    // test 3
    test = {
        "ACGTACGTACGTACGTACGTACGTACGTACGAAGCACGTACGTACGTACGTACGTACGTACGTACGT",
        "CGGC",
        {40,40,40,40},
        {90,90,90,90,90,90,90,90,
         90,90,90,90,90,90,90,90,
         90,90,90,90,90,90,90,90,
         90,90,90,90,90,90,90,70,
         90,90,90,90,90,90,90,90,
         90,90,90,90,90,90,90,90,
         90,90,90,90,90,90,90,90,
         90,90,90,90,90,90,90,90,
         90,90,90},
        1,
        4
    };
    expected_alignment = {
        71,
        29,
        "CGAAGC",
        "CG--GC"
    };
    CHECK_TEST(test, sse_short_hmm)
    CHECK_ALIGNER(test, sse_short_hmm, expected_alignment)
    CHECK_TEST(test, sse_int_hmm)
    CHECK_ALIGNER(test, sse_int_hmm, expected_alignment)
    
    // test 4
    test = {
        "CGAAGCACGTACGTACGTAACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGTACGT",
        "CGGC",
        {40,40,40,40},
        {90,90,70,90,90,90,90,90,
         90,90,90,90,90,90,90,90,
         90,90,90,90,90,90,90,90,
         90,90,90,90,90,90,90,90,
         90,90,90,90,90,90,90,90,
         90,90,90,90,90,90,90,90,
         90,90,90,90,90,90,90,90,
         90,90,90,90,90,90,90,90,
         90,90,90},
        1,
        4
    };
    expected_alignment = {
        71,
        0,
        "CGAAGC",
        "CG--GC"
    };
    CHECK_TEST(test, sse_short_hmm)
    CHECK_ALIGNER(test, sse_short_hmm, expected_alignment)
    CHECK_TEST(test, sse_int_hmm)
    CHECK_ALIGNER(test, sse_int_hmm, expected_alignment)
    
    // test 5
    test = {
        "CCCCACGTCCCCACGTCCCCACGTCCCCACGTATATATATATATATGGGGACGTGGGGACGTGGGGACGTGGGGACGT",
        "CCCCACGTGGGACGT",
        {40,40,40,40,40,40,40,40,40,40,40,40,40,40,40},
        {90,90,90,90,90,90,90,90,
         90,90,90,90,90,90,90,90,
         90,90,90,90,90,90,90,90,
         90,90,90,90,90,90,90,90,
         70,90,90,90,90,90,90,90,
         90,90,90,90,90,90,90,90,
         90,90,90,90,90,90,90,90,
         90,90,90,90,90,90,90,90,
         90,90,90,90,90,90,90,90,
         90,90,90,90,90,90},
        1,
        4
    };
    expected_alignment = {
        84,
        24,
        "CCCCACGTATATATATATATATGGGGACGT",
        "CCCCACGT---------------GGGACGT"
    };
    CHECK_TEST(test, sse_short_hmm)
    CHECK_ALIGNER(test, sse_short_hmm, expected_alignment)
    CHECK_TEST(test, sse_int_hmm)
    CHECK_ALIGNER(test, sse_int_hmm, expected_alignment)
}

#ifdef __AVX2__
BOOST_AUTO_TEST_CASE(avx2_check_alignments)
{
    AVX2PairHMM<16, short> hmm;
    TestCase test;
    Alignment expected_alignment;
    
    // test 1
    test = {
        "ACGTACGTACGTACGTACGTACGTACGTACGAAAA",
        "AAAA",
        {40,40,40,40},
        {10,10,10,10,10,10,10,10,
         10,10,10,10,10,10,10,10,
         10,10,10,10,10,10,10,10,
         10,10,10,10,10,10,10,10,
         10,10,10},
        1,
        4
    };
    expected_alignment = {
        0,
        31,
        "AAAA",
        "AAAA"
    };
    CHECK_TEST(test, hmm)
    CHECK_ALIGNER(test, hmm, expected_alignment);
    
    // test 2
    test = {
        "ACGTACGTACGTACGTACGTACGTACGTACGAATA",
        "AAAA",
        {40,40,40,40},
        {90,90,90,90,90,90,90,90,
         90,90,90,90,90,90,90,90,
         90,90,90,90,90,90,90,90,
         90,90,90,90,90,90,90,90,
         90,90,90},
        1,
        4
    };
    expected_alignment = {
        40,
        31,
        "AATA",
        "AAAA"
    };
    CHECK_TEST(test, hmm)
    CHECK_ALIGNER(test, hmm, expected_alignment);
    
    // test 3
    test = {
        "ACGTACGTACGTACGAAGCACGTACGTACGTACGT",
        "CGGC",
        {40,40,40,40},
        {90,90,90,90,90,90,90,90,
         90,90,90,90,90,90,90,70,
         90,90,90,90,90,90,90,90,
         90,90,90,90,90,90,90,90,
         90,90,90},
        1,
        4
    };
    expected_alignment = {
        71,
        13,
        "CGAAGC",
        "CG--GC"
    };
    CHECK_TEST(test, hmm)
    CHECK_ALIGNER(test, hmm, expected_alignment);
    
    // test 4
    test = {
        "CGAAGCACGTACGTACGTAACGTACGTACGTACGT",
        "CGGC",
        {40,40,40,40},
        {90,90,70,90,90,90,90,90,
         90,90,90,90,90,90,90,90,
         90,90,90,90,90,90,90,90,
         90,90,90,90,90,90,90,90,
         90,90,90},
        1,
        4
    };
    expected_alignment = {
        71,
        0,
        "CGAAGC",
        "CG--GC"
    };
    CHECK_TEST(test, hmm)
    CHECK_ALIGNER(test, hmm, expected_alignment);

    // test 5
    test = {
        "CCCCACGTCCCCACGTATATATATATATATGGGGACGTGGGGACGT",
        "CCCCACGTGGGACGT",
        {40,40,40,40,40,40,40,40,40,40,40,40,40,40,40},
        {90,90,90,90,90,90,90,90,
         90,90,90,90,90,90,90,90,
         70,90,90,90,90,90,90,90,
         90,90,90,90,90,90,90,90,
         90,90,90,90,90,90,90,90,
         90,90,90,90,90,90},
        1,
        4
    };
    expected_alignment = {
        84,
        8,
        "CCCCACGTATATATATATATATGGGGACGT",
        "CCCCACGT---------------GGGACGT"
    };
    CHECK_TEST(test, hmm)
    CHECK_ALIGNER(test, hmm, expected_alignment);
}
#endif /* __AVX2__ */


// Speed tests


const int iterations = 1000000;


BOOST_AUTO_TEST_CASE(sse2_short_band_size_8_speed)
{
    SSE2PairHMM<8, short> sse_short_hmm;

    CHECK_TEST(band8_speed_test, sse_short_hmm)
    CHECK_ALIGNER(band8_speed_test, sse_short_hmm, band8_speed_expected_alignment)
    CHECK_SPEED(band8_speed_test, sse_short_hmm, band8_speed_expected_alignment, iterations)
};


BOOST_AUTO_TEST_CASE(sse2_int_band_size_8_speed)
{
    SSE2PairHMM<8, int> sse_int_hmm;

    CHECK_TEST(band8_speed_test, sse_int_hmm)
    CHECK_ALIGNER(band8_speed_test, sse_int_hmm, band8_speed_expected_alignment)
    CHECK_SPEED(band8_speed_test, sse_int_hmm, band8_speed_expected_alignment, iterations)
};


BOOST_AUTO_TEST_CASE(sse2_short_band_size_16_speed)
{
    SSE2PairHMM<16, short> sse_short_hmm;
    CHECK_TEST(band16_speed_test, sse_short_hmm)
    CHECK_ALIGNER(band16_speed_test, sse_short_hmm, band8_speed_expected_alignment)
    CHECK_SPEED(band16_speed_test, sse_short_hmm, band8_speed_expected_alignment, iterations)
};


BOOST_AUTO_TEST_CASE(sse2_int_band_size_16_speed)
{
    SSE2PairHMM<16, int> sse_int_hmm;
    CHECK_TEST(band16_speed_test, sse_int_hmm)
    CHECK_ALIGNER(band16_speed_test, sse_int_hmm, band8_speed_expected_alignment)
    CHECK_SPEED(band16_speed_test, sse_int_hmm, band8_speed_expected_alignment, iterations)
};


BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

} // namespace test
} // namespace octopus
