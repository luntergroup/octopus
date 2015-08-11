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
#include "htslib_bcf_facade.h"
#include "vcf_header.h"
#include "vcf_record.h"
#include "vcf_reader.h"

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
    auto v6 = header.get_typed_value("INFO", "SOMATIC", "0");
    
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
    
    //cout << header << endl;
    
    GenomicRegion region {"X", 2000000, 2001000};
    
    auto start = std::chrono::system_clock::now();
    
    auto records = reader.fetch_records();
    
    auto end = std::chrono::system_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start).count();
    
    cout << duration << endl;
    
    cout << records.size() << endl;
    
//    for (const auto& record : records) {
//        auto dp = get_typed_info_values("DP", record.get_info_value("DP"), header);
//        cout << dp.front() << endl;
//    }
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
