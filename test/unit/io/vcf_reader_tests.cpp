// Copyright (c) 2016 Daniel Cooke
// Use of this source code is governed by the MIT license that can be found in the LICENSE file.

#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include <iostream>
#include <string>

#include "test_common.hpp"
#include <basics/genomic_region.hpp>
#include <core/types/variant.hpp>
#include "vcf.hpp"

using std::cout;
using std::endl;

namespace octopus { namespace test {

BOOST_AUTO_TEST_SUITE(Components)
BOOST_AUTO_TEST_SUITE(IO)

BOOST_AUTO_TEST_CASE(can_read_vcf_files)
{
    BOOST_REQUIRE(test_file_exists(sample_vcf));
    
    VcfReader reader {sample_vcf};
    
    BOOST_REQUIRE(reader.is_open());
    
    const GenomicRegion region {"X", 2000000, 2001000};
    
    auto records = reader.fetch_records(region);
    
    //for (const auto& record : records) cout << record << endl;
}

BOOST_AUTO_TEST_CASE(can_use_vcf_types)
{
    BOOST_REQUIRE(test_file_exists(sample_vcf));
    
    VcfReader reader {sample_vcf};
    
    BOOST_REQUIRE(reader.is_open());
    
    const auto header = reader.fetch_header();
    
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
    
    BOOST_REQUIRE(writer.is_open());
    
    const auto header = VcfHeader::Builder().add_contig("TEST").build_once();
    const auto record = VcfRecord::Builder().set_chrom("TEST").set_id("TEST").set_pos(0).set_qual(60).set_ref("A").set_alt("C").build_once();
    
    writer.write(header);
    
    BOOST_REQUIRE(writer.is_open());
    
    writer.write(record);
    writer.close();
    
    BOOST_CHECK(test_file_exists(test_out_vcf));
    
    remove_test_file(test_out_vcf);
    
    BOOST_CHECK(!test_file_exists(test_out_vcf));
}

BOOST_AUTO_TEST_CASE(can_write_vcfgz_files)
{
    VcfWriter writer {test_out_vcfgz};
    
    BOOST_REQUIRE(writer.is_open());
    
    const auto header = VcfHeader::Builder().add_contig("TEST").build_once();
    const auto record = VcfRecord::Builder().set_chrom("TEST").set_id("TEST").set_pos(0).set_qual(60).set_ref("A").set_alt("C").build_once();
    
    writer.write(header);
    
    BOOST_REQUIRE(writer.is_open());
    
    writer.write(record);
    writer.close();
    
    BOOST_CHECK(test_file_exists(test_out_vcfgz));
    
    remove_test_file(test_out_vcfgz);
    
    BOOST_CHECK(!test_file_exists(test_out_vcfgz));
}

BOOST_AUTO_TEST_CASE(can_write_bcf_files)
{
    VcfWriter writer {test_out_bcf};
    
    BOOST_REQUIRE(writer.is_open());
    
    const auto header = VcfHeader::Builder().add_contig("TEST").build_once();
    const auto record = VcfRecord::Builder().set_chrom("TEST").set_id("TEST").set_pos(0).set_qual(60).set_ref("A").set_alt("C").build_once();
    
    writer.write(header);
    
    BOOST_REQUIRE(writer.is_open());
    
    writer.write(record);
    writer.close();
    
    BOOST_CHECK(test_file_exists(test_out_bcf));
    
    remove_test_file(test_out_bcf);
    
    BOOST_CHECK(!test_file_exists(test_out_bcf));
}

// This test seems to screw up error reporting for the others so will need to test independently
BOOST_AUTO_TEST_CASE(can_write_vcf_to_stdout)
{
//    VcfWriter writer {"-"};
//    BOOST_CHECK(writer.is_open()); // the only practical thing we can check here
}

BOOST_AUTO_TEST_CASE(vcfgz_files_can_be_indexed)
{
    
}

BOOST_AUTO_TEST_CASE(bcf_files_can_be_indexed)
{
    
}

BOOST_AUTO_TEST_SUITE_END() // IO
BOOST_AUTO_TEST_SUITE_END() // Components

} // namespace octopus
} // namespace test
