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
#include <algorithm>

#include "test_common.h"
#include "genomic_region.h"
#include "variant.h"
#include "variant_file_reader.h"
#include "variant_file_factory.h"
#include "vcf_header.h"
#include "vcf_type.h"
#include "vcf_record.h"
#include "vcf_reader.h"
#include "vcf_writer.h"
#include "vcf_parser.h"

using std::cout;
using std::endl;

std::vector<GenomicRegion> get_batch_regions(const GenomicRegion& region, const VcfReader& reader, std::size_t max_batch_size)
{
    std::vector<GenomicRegion> result {};
    
    if (reader.num_records(region) > max_batch_size) {
        // umm?
    } else {
        result.push_back(region);
    }
    
    return result;
}

std::vector<Variant> fetch_variants(const GenomicRegion& region, VcfReader& reader)
{
    std::vector<Variant> result {};
    result.reserve(reader.num_records(region));
    
    std::size_t max_batch_size {10000};
    
    auto batches = get_batch_regions(region, reader, max_batch_size);
    
    for (const auto& batch : batches) {
        auto vcf_records = reader.fetch_records(batch, VcfReader::Unpack::AllButSamples);
        for (const auto& record : vcf_records) {
            for (const auto& alt_allele : record.get_alt_alleles()) {
                result.emplace_back(record.get_chromosome_name(), record.get_position(), record.get_ref_allele(), alt_allele);
            }
        }
    }
    
    return result;
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
}

void test3()
{
    GenomicRegion region {"X", 2000000, 2010000};
    
    VcfReader reader {sample_vcf};
    VcfParser parser {"/Users/dcooke/Genomics/sample_vcf/CEU.low_coverage.2010_07.xchr.genotypes.vcf"};
    
    //auto header = reader.fetch_header();
    auto header = VcfHeader::Builder(parser.fetch_header()).add_contig("X").build_once();
    //auto header = parser.fetch_header();
    
    //auto records = reader.fetch_records(region);
    auto records = parser.fetch_records(region);
    
    VcfWriter writer {"/Users/dcooke/test.vcf.gz"};
    writer.write(header);
    for (const auto& record : records) {
        writer.write(record);
    }
    
    //    auto variants = fetch_variants(region, reader);
    //
    //    for (const auto& variant : variants) {
    //        cout << variant << endl;
    //    }
    
    //    cout << reader.num_records(region) << endl;
    //    cout << parser.num_records(region) << endl;
    
    //auto records = reader.fetch_records(VcfReader::Unpack::AllButSamples);
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
