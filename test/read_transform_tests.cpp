//
//  read_transform_tests.cpp
//  Octopus
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
#include "genomic_region.hpp"
#include "read_manager.hpp"
#include "read_transform.hpp"
#include "read_transformations.hpp"

BOOST_AUTO_TEST_SUITE(Components)

BOOST_AUTO_TEST_CASE(read_transform_test)
{
    BOOST_REQUIRE(test_file_exists(NA12878_low_coverage));
    
    ReadManager read_manager {NA12878_low_coverage};
    
    auto sample = read_manager.get_samples().at(0);
    
    GenomicRegion region1 {"Y", 3000000, 3000010};
    
    auto reads = read_manager.fetch_reads(sample, region1);
    
    BOOST_REQUIRE(std::is_sorted(std::cbegin(reads), std::cend(reads)));
    
    BOOST_REQUIRE(reads.size() == 4);
    
    const auto& a_read = reads[0];
    
    BOOST_CHECK(is_back_soft_clipped(a_read.get_cigar_string()));
    
    Octopus::ReadTransform transformer {};
    transformer.register_transform(Octopus::ReadTransforms::trim_adapters());
    transformer.register_transform(Octopus::ReadTransforms::trim_soft_clipped());
    
    transform_reads(reads, transformer);
    
    BOOST_CHECK(std::all_of(a_read.get_qualities().rbegin(), a_read.get_qualities().rbegin() + 13,
                            [] (const auto q) { return q == 0; }));
    
    GenomicRegion region2 {"18", 389260, 389361};
    
    reads = read_manager.fetch_reads(sample, region2);
    
    transform_reads(reads, transformer);
    
    // TODO
}

BOOST_AUTO_TEST_SUITE_END()
