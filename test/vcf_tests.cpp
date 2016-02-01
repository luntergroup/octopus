//
//  vcf_tests.cpp
//  Octopus
//
//  Created by Daniel Cooke on 15/05/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include <iostream>
#include <string>

#include "test_common.hpp"
#include "genomic_region.hpp"
#include "variant.hpp"
#include "vcf.hpp"

using std::cout;
using std::endl;

BOOST_AUTO_TEST_CASE(can_read_vcf_files)
{
    BOOST_REQUIRE(test_file_exists(sample_vcf));
    
    VcfReader vcf_reader {sample_vcf};
    
    GenomicRegion region {"X", 2000000, 2001000};
    
    auto records = vcf_reader.fetch_records(region);
    
    //for (const auto& record : records) cout << record << endl;
}

BOOST_AUTO_TEST_CASE(can_use_vcf_types)
{
    BOOST_REQUIRE(test_file_exists(sample_vcf));
    
    VcfReader reader {sample_vcf};
    
    auto header = reader.fetch_header();
    
    auto v1 = get_typed_info_value(header, "AA", "ACGT");
    auto v2 = get_typed_info_value(header, "DP", "110");
    auto v3 = get_typed_info_value(header, "AF", "0.9");
    //auto v4 = get_typed_info_value(header, "SOMATIC", "1");
    auto v5 = get_typed_info_value(header, "AA", "TGCA");
    //auto v6 = get_typed_info_value(header, "SOMATIC", "0");
}

BOOST_AUTO_TEST_CASE(can_write_vcf_files)
{
    VcfWriter writer {test_out_vcf};
    
    auto header  = VcfHeader::Builder().add_basic_field("TEST", "TEST").build_once();
    auto record  = VcfRecord::Builder().set_chromosome("TEST").set_id("TEST").set_position(0).set_quality(60).set_ref_allele("A").set_alt_allele("C").build_once();
    
    writer.write(header);
    writer.write(record);
    
    remove_test_file(test_out_vcf);
}

BOOST_AUTO_TEST_CASE(can_write_vcf_to_stdout)
{
    // TODO
//    VcfWriter writer {"-"};
//    
//    auto header  = VcfHeader::Builder().add_basic_field("TEST", "TEST").build_once();
//    auto record  = VcfRecord::Builder().set_chromosome("TEST").set_id("TEST").set_position(0).set_quality(60).set_ref_allele("A").set_alt_allele("C").build_once();
//    
//    writer.write(header);
//    writer.write(record);
}
