// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#include <boost/test/unit_test.hpp>

#include <stdexcept>
#include <vector>
#include <iterator>
#include <numeric>
#include <algorithm>
#include <cstdint>
#include <cmath>
#include <functional>

#include <boost/test/floating_point_comparison.hpp>

#include <basics/phred.hpp>

namespace octopus { namespace test {

namespace { constexpr double tolerance {0.000001}; }

BOOST_AUTO_TEST_SUITE(basics)
BOOST_AUTO_TEST_SUITE(phred)

BOOST_AUTO_TEST_CASE(phreds_must_be_positive)
{
    BOOST_CHECK_THROW((Phred<> {-1}), std::exception);
    BOOST_CHECK_NO_THROW((Phred<> {0}));
    BOOST_CHECK_NO_THROW((Phred<> {-0}));
}

BOOST_AUTO_TEST_CASE(phred_are_the_size_of_their_underlying_type)
{
    BOOST_CHECK_EQUAL(sizeof(Phred<signed char> {}),   sizeof(signed char));
    BOOST_CHECK_EQUAL(sizeof(Phred<unsigned char> {}), sizeof(unsigned char));
    BOOST_CHECK_EQUAL(sizeof(Phred<int> {}),           sizeof(int));
    BOOST_CHECK_EQUAL(sizeof(Phred<unsigned> {}),      sizeof(unsigned));
    BOOST_CHECK_EQUAL(sizeof(Phred<float> {}),         sizeof(float));
    BOOST_CHECK_EQUAL(sizeof(Phred<double> {}),        sizeof(double));
    BOOST_CHECK_EQUAL(sizeof(Phred<std::uint8_t> {}),  sizeof(std::uint8_t));
    BOOST_CHECK_EQUAL(sizeof(Phred<std::int8_t> {}),   sizeof(std::int8_t));
    BOOST_CHECK_EQUAL(sizeof(Phred<std::uint16_t> {}), sizeof(std::uint16_t));
    BOOST_CHECK_EQUAL(sizeof(Phred<std::int16_t> {}),  sizeof(std::int16_t));
}

BOOST_AUTO_TEST_CASE(phreds_can_be_converted_to_probabilities)
{
    std::vector<double> v(100);
    
    std::iota(std::begin(v), std::end(v), 0.0);
    
    for (auto x : v) {
        Phred<> phred {x};
        BOOST_CHECK_EQUAL(phred.score(), x);
        //BOOST_CHECK_CLOSE(phred.probability_false(), std::pow(10.0, -x / 10.0), tolerance);
    }
}

BOOST_AUTO_TEST_CASE(phreds_can_be_constructed_with_probabilities)
{
    std::vector<double> v(100, 0.1);
    
    v.front() = 1;
    std::partial_sum(std::next(std::cbegin(v)), std::cend(v), std::next(std::begin(v)), std::multiplies<> {});
    
    double c {0};
    for (auto p : v) {
        auto phred = probability_to_phred(p);
        BOOST_CHECK_CLOSE(phred.score(), c++, tolerance);
        //BOOST_CHECK_CLOSE(phred.probability_false(), p, tolerance);
    }
}

BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()

} // namespace test
} // namespace octopus
