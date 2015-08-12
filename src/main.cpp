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
#include <chrono>

#include "test_common.h"
#include "genomic_region.h"
#include "variant.h"
#include "variant_file_reader.h"
#include "variant_file_factory.h"
#include "vcf_header.h"
#include "vcf_record.h"
#include "vcf_reader.h"
#include "vcf_writer.h"

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
    VcfReader reader {sample_vcf};
    
    auto header = reader.fetch_header();
    
    auto v1 = get_typed_info_value(header, "AA", "ACGT");
    auto v2 = get_typed_info_value(header, "DP", "110");
    auto v3 = get_typed_info_value(header, "AF", "0.9");
    auto v4 = get_typed_info_value(header, "SOMATIC", "1");
    auto v5 = get_typed_info_value(header, "AA", "TGCA");
    auto v6 = get_typed_info_value(header, "SOMATIC", "0");
    
    //v3 += 5;
    
//    cout << v2 << endl;
//    cout << v3 << endl;
//    cout << (v2 + v3) << endl;
    
//    auto x = static_cast<double>(v2);
//    auto y = static_cast<double>(v3);
//    
//    mu::Parser parser {};
//    
//    parser.DefineVar("x", &x);
//    parser.DefineVar("y", &y);
//    
//    parser.SetExpr("x + y");
//    
//    cout << parser.Eval() << endl;
}

void test3()
{
    VcfReader reader {sample_vcf};
    
    auto header = reader.fetch_header();
    
    GenomicRegion region {"X", 2000000, 2001000};
    
    auto records = reader.fetch_records(region);
    
    VcfWriter writer {"/Users/danielcooke/test.vcf"};
    
    writer.write(header);
    //writer.write(records.front());
}

int main(int argc, const char **argv)
{
    test3();
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
