//
//  mappable_set_tests.cpp
//  Octopus
//
//  Created by Daniel Cooke on 19/07/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include <iostream>
#include <chrono>
#include <algorithm>

#include "test_common.hpp"
#include "genomic_region.hpp"
#include "read_manager.hpp"
#include "mappable_set.hpp"
#include "read_utils.hpp"
#include "read_filters.hpp"
#include "context_iterators.hpp"
#include "mappable_map.hpp"

using std::cout;
using std::endl;

BOOST_AUTO_TEST_SUITE(Components)

BOOST_AUTO_TEST_CASE(MappableSet_works_like_vector)
{
    
}

BOOST_AUTO_TEST_SUITE_END()
