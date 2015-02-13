//
//  read_reader_test.cpp
//  Octopus
//
//  Created by Daniel Cooke on 11/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "catch.hpp"

#include <iostream>
#include <string>

#include "htslib_facade.h"

TEST_CASE("read_reader_open_test", "[read_reader]")
{
    std::string homedir {getenv("HOME")};
    std::string the_bam_name {"HG00101.mapped.ILLUMINA.bwa.GBR.low_coverage.20130415.bam"};
    std::string bam_file {homedir + "/Genomics/Illumina/" + the_bam_name};
    
    HtslibFacade a_reader {bam_file};
    auto reads = a_reader.fetch_reads(GenomicRegion("10", 1000000, 1000100));
    
    for (auto read : reads) {
        std::cout << read << std::endl;
    }
    
    REQUIRE(true);
}
