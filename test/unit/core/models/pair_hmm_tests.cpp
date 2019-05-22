// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include <boost/test/unit_test.hpp>

#include <string>
#include <vector>
#include <cstddef>
#include <algorithm>

#include "core/models/pairhmm/simd_pair_hmm.hpp"

namespace octopus { namespace test {

using namespace octopus::hmm::simd;

BOOST_AUTO_TEST_SUITE(core)
BOOST_AUTO_TEST_SUITE(model)

auto add_pad(std::string seq, const int pad = 8)
{
    return std::string(pad, 'N') + std::move(seq) + std::string(pad, 'N');
}

auto add_pad(std::vector<std::int8_t> quals, const int pad = 8, std::int8_t val = 50)
{
    quals.resize(quals.size() + 2 * pad, 50);
    std::rotate(quals.begin(), quals.begin() + pad, quals.end());
    return quals;
}

auto align_helper(const std::string& target,
                  const std::string& query,
                  const std::vector<std::int8_t>& qualities,
                  int gap_open, int gap_extend, int nuc_prior,
                  int pad = 8)
{
    const auto padded_target = add_pad(target, pad);
    const auto padded_qualities = add_pad(qualities, pad);
    return align(padded_target.data(), query.data(), padded_qualities.data(),
                 static_cast<int>(padded_target.size()), static_cast<int>(query.size()),
                 gap_open, gap_extend, nuc_prior);
}

BOOST_AUTO_TEST_CASE(correct_alignment_score)
{
    std::string target, query;
    std::vector<std::int8_t> qualities;
    int score;
    
    target = "AAAAAAAAAAAAAAAAAA";
    query = "AAAAAAAAAAAAAAAAAA";
    qualities.assign(query.size(), 30);
    score = align_helper(target, query, qualities, 45, 3, 50);
    BOOST_CHECK_EQUAL(score, 0);
    
    target = "AAAAAAACCCCAAAAAAA";
    query = "AAAAAAACCCCAAAAAAA";
    qualities.assign(query.size(), 30);
    score = align_helper(target, query, qualities, 45, 3, 50);
    BOOST_CHECK_EQUAL(score, 0);
    
    target = "AAAAAAACCCCAAAAAAA";
    query = "AAAAAAACACCAAAAAAA";
    qualities.assign(query.size(), 30);
    score = align_helper(target, query, qualities, 45, 3, 50);
    BOOST_CHECK_EQUAL(score, 30);
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

} // namespace test
} // namespace octopus
