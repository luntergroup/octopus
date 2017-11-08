// Copyright (c) 2017 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.
 
#include <stdio.h>
#include <boost/test/unit_test.hpp>
#include <boost/test/floating_point_comparison.hpp>
 
#include "utils/maths.hpp"
 
namespace octopus { namespace test {
 
static constexpr double tolerance {1e-10};
 
BOOST_AUTO_TEST_SUITE(utils)
BOOST_AUTO_TEST_SUITE(maths)
 
BOOST_AUTO_TEST_CASE(log_sum_exp_handles_edge_cases)
{
    using octopus::maths::log_sum_exp;
    static constexpr double zero {0.0};
    static constexpr double lnHalf {-0.6931471805599453};
    BOOST_CHECK_CLOSE(log_sum_exp(lnHalf, lnHalf), zero, tolerance);
    BOOST_CHECK_CLOSE(log_sum_exp(zero, zero), -lnHalf, tolerance);
}
 
BOOST_AUTO_TEST_SUITE_END()
BOOST_AUTO_TEST_SUITE_END()
 
} // namespace test
} // namespace octopus
