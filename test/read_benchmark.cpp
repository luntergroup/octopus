//
//  read_benchmark.cpp
//  Octopus
//
//  Created by Daniel Cooke on 14/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include <iostream>

#include "catch.hpp"
#include "test_common.h"
#include "benchmark_utils.h"

#include "htslib_facade.h"
#include "read_factory.h"

TEST_CASE("read_benchmark", "[benchmark]")
{
    HtslibFacade a_reader {human_1000g_bam};
    
    auto f_read = [&a_reader] () {
        a_reader.fetch_reads(GenomicRegion("10", 1000000, 1010000));
    };
    
    auto without_vptr = benchmark<std::chrono::microseconds>(f_read, 10).count();
    
    ReadFactory a_read_factory(std::vector<std::string> {human_1000g_bam});
    
    auto sample_ids = a_read_factory.get_sample_ids();
    
    auto the_sample_id = sample_ids.at(0);
    
    auto f_factory = [&a_read_factory, &the_sample_id] () {
        a_read_factory.fetch_reads(the_sample_id, GenomicRegion("10", 1000000, 1010000));
    };
    
    auto with_vptr = benchmark<std::chrono::microseconds>(f_factory, 10).count();
    
    //std::cout << "Without vptr: " << without_vptr << "us" << std::endl;
    //std::cout << "With vptr: " << with_vptr << "us" << std::endl;
    
    REQUIRE(true);
}
