//
//  read_transform_tests.cpp
//  octopus
//
//  Created by Daniel Cooke on 08/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include <iostream>
#include <vector>
#include <algorithm>
#include <iterator>

#include "test_common.hpp"
#include <basics/genomic_region.hpp>
#include <io/read/read_manager.hpp>
#include "read_transformer.hpp"
#include "read_transform.hpp"

using std::cout;
using std::endl;

BOOST_AUTO_TEST_SUITE(Components)

//BOOST_AUTO_TEST_CASE(read_transform_test)
//{
//    BOOST_REQUIRE(test_file_exists(NA12878_low_coverage));
//    
//    ReadManager read_manager {NA12878_low_coverage};
//    
//    const auto sample = read_manager.samples().front();
//    
//    GenomicRegion region1 {"10", 1'000'000, 1'000'100};
//    
//    auto reads = read_manager.fetch_reads(sample, region1);
//    
//    BOOST_REQUIRE(std::is_sorted(std::cbegin(reads), std::cend(reads)));
//    BOOST_REQUIRE(reads.size() == 7);
//    
//    const auto& a_read = reads[0];
//    
//    BOOST_CHECK(is_back_soft_clipped(a_read.cigar_string()));
//    
//    octopus::ReadTransform transformer {};
//    transformer.register_transform(octopus::ReadTransforms::trim_adapters());
//    transformer.register_transform(octopus::ReadTransforms::trim_soft_clipped());
//    
//    transform_reads(reads, transformer);
//    
//    BOOST_CHECK(std::all_of(a_read.qualities().rbegin(), a_read.qualities().rbegin() + 13,
//                            [] (const auto q) { return q == 0; }));
//    
//    GenomicRegion region2 {"3", 100'000, 100'100};
//    
//    reads = read_manager.fetch_reads(sample, region2);
//    
//    BOOST_REQUIRE(std::is_sorted(std::cbegin(reads), std::cend(reads)));
//    BOOST_REQUIRE(reads.size() == 21);
//    
//    transform_reads(reads, transformer);
//    
//    // TODO
//}

BOOST_AUTO_TEST_SUITE_END()
