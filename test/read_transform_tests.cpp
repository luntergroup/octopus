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
#include <algorithm> // std::sort
#include <iterator>  // std::back_inserter

#include "test_common.hpp"
#include "genomic_region.hpp"
#include "read_manager.hpp"
#include "read_transform.hpp"
#include "read_transformations.hpp"

BOOST_AUTO_TEST_SUITE(Components)

BOOST_AUTO_TEST_CASE(read_transform_test)
{
    ReadManager a_read_manager {HG00101};
    
    auto sample_ids = a_read_manager.get_samples();
    
    auto the_sample_id = sample_ids.at(0);
    
    GenomicRegion a_region {"Y", 3000000, 3000010};
    
    auto some_reads = a_read_manager.fetch_reads(the_sample_id, a_region);
    
    if (!std::is_sorted(some_reads.begin(), some_reads.end())) {
        std::sort(some_reads.begin(), some_reads.end());
    }
    
    BOOST_CHECK(some_reads.size() == 4);
    
    const auto& a_read = some_reads[0];
    
    BOOST_CHECK(is_back_soft_clipped(a_read.get_cigar_string()));
    
    Octopus::ReadTransform transformer {};
    transformer.register_transform(Octopus::ReadTransforms::trim_adapters());
    transformer.register_transform(Octopus::ReadTransforms::trim_soft_clipped());
    transformer.transform_reads(some_reads.begin(), some_reads.end());
    
    BOOST_CHECK(std::all_of(a_read.get_qualities().rbegin(), a_read.get_qualities().rbegin() + 13,
                            [] (auto q) { return q == 0; }));
    
    GenomicRegion another_region {"18", 389260, 389361};
    
    some_reads = a_read_manager.fetch_reads(the_sample_id, another_region);
    
    //auto& read_with_adapter = some_reads[6];
    
    //std::cout << read_with_adapter << std::endl;
    
    transformer.transform_reads(some_reads.begin(), some_reads.end());
    
}

BOOST_AUTO_TEST_SUITE_END()
