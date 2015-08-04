//
//  main.cpp
//  Octopus
//
//  Created by Daniel Cooke on 03/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

//#define BOOST_TEST_DYN_LINK
//#define BOOST_TEST_MODULE Main
//#include <boost/test/unit_test.hpp>

//#include <iostream>
//#include <cstdlib>
//
//#include "program_options.h"
//#include "octopus.h"

#include <iostream>
#include <string>

#include "test_common.h"
#include "genomic_region.h"
#include "variant.h"
#include "variant_file_reader.h"
#include "variant_file_factory.h"
#include "htslib_bcf_facade.h"
#include "vcf_record.h"

using std::cout;
using std::endl;

void test()
{
    HtslibBcfFacade vcf_reader {sample_vcf};
    
    GenomicRegion region {"X", 2000000, 2001000};
    
    auto records = vcf_reader.fetch_records(region);
    
    for (const auto& record : records) cout << record << endl;
}

int main(int argc, const char **argv)
{
    test();
//    auto options = Octopus::parse_options(argc, argv);
//    
//    if (options.second) {
//        Octopus::run_octopus(options.first);
//        std::cout << "finished running Octopus" << std::endl;
//    } else {
//        std::cout << "did not run Octopus" << std::endl;
//    }
    
    return EXIT_SUCCESS;
}
