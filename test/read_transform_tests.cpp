//
//  read_transform_tests.cpp
//  Octopus
//
//  Created by Daniel Cooke on 08/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "catch.hpp"

#include <iostream>
#include <vector>
#include <algorithm> // std::sort
#include <iterator>  // std::back_inserter

#include "test_common.h"
#include "genomic_region.h"
#include "read_manager.h"
#include "read_transform.h"
#include "read_transformations.h"

TEST_CASE("read_transform_test", "[read_transform]")
{
    
    ReadManager a_read_manager(std::vector<std::string> {human_1000g_bam});
    
    auto sample_ids = a_read_manager.get_sample_ids();
    
    auto the_sample_id = sample_ids.at(0);
    
    GenomicRegion a_region {"Y", 3000000, 3000010};
    
    std::vector<AlignedRead> some_reads = a_read_manager.fetch_reads(the_sample_id, a_region);
    
    if (!std::is_sorted(some_reads.begin(), some_reads.end())) {
        std::sort(some_reads.begin(), some_reads.end());
    }
    
    REQUIRE(some_reads.size() == 4);
    
    const auto& a_read = some_reads[0];
    
    REQUIRE(is_back_soft_clipped(a_read.get_cigar_string()));
    
    ReadTransform a_read_transform {};
    a_read_transform.register_transform(trim_adapters);
    a_read_transform.register_transform(trim_soft_clipped);
    a_read_transform.transform_reads(some_reads.begin(), some_reads.end());
    
    REQUIRE(std::all_of(a_read.get_qualities().rbegin(), a_read.get_qualities().rbegin() + 13,
                        [] (auto q) { return q == 0; }));
    
    GenomicRegion another_region {"18", 389260, 389361};
    
    some_reads = a_read_manager.fetch_reads(the_sample_id, another_region);
    
    auto& read_with_adapter = some_reads[6];
    
    //std::cout << read_with_adapter << std::endl;
    
    a_read_transform.transform_reads(some_reads.begin(), some_reads.end());
    
    
}
