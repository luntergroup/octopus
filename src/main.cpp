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
#include "vcf_header.h"
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

void test2()
{
    VcfHeader header {"##INFO=<ID=AA,Type=String>\n##INFO=<ID=DP,Type=Integer>\n##INFO=<ID=SOMATIC,Type=Flag>\n##INFO=<ID=AF,Type=Float>"};
    
    cout << header << endl;
    
    auto v1 = header.get_typed_value("INFO", "AA", "ACGT");
    auto v2 = header.get_typed_value("INFO", "DP", "110");
    auto v3 = header.get_typed_value("INFO", "AF", "0.9");
    auto v4 = header.get_typed_value("INFO", "SOMATIC", "1");
    auto v5 = header.get_typed_value("INFO", "AA", "TGCA");
    
//    cout << v1 << endl;
//    cout << v2 << endl;
//    cout << v3 << endl;
//    cout << v4 << endl;
    
    auto v6 = v2 * v3;
    cout << v6 << endl;
    
    auto v7 = v1 + v5;
    cout << v7 << endl;
    
    cout << (v2 == v3) << endl;
    cout << (v2 < v3) << endl;
    cout << (v2 > v3) << endl;
    cout << (v2 <= v3) << endl;
    cout << (v2 >= v3) << endl;
    cout << (v6 > v2) << endl;
}

int main(int argc, const char **argv)
{
    test2();
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
