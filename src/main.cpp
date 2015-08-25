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
#include "reference_genome.h"
#include "vcf_header.h"
#include "vcf_type.h"
#include "vcf_record.h"
#include "vcf_reader.h"
#include "vcf_writer.h"
#include "vcf_parser.h"
#include "read_manager.h"
#include "mappable_set.h"
#include "candidate_variant_generator.h"
#include "alignment_candidate_variant_generator.h"

#include "sequence_utils.h"

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
    
    //VcfReader reader {sample_vcf};
    VcfReader reader {"/Users/danielcooke/Genomics/sample_vcf/CEU.low_coverage.2010_07.xchr.genotypes.vcf"};
    
    auto header = reader.fetch_header();
    //auto header = VcfHeader::Builder(reader.fetch_header()).add_contig("X").build_once();
    
    auto records = reader.fetch_records(region);
    
    cout << records.front() << endl;
    
//    VcfWriter writer {"/Users/danielcooke/test.vcf.gz"};
//    writer.write(header);
//    for (const auto& record : records) {
//        writer.write(record);
//    }
}

void test4()
{
    GenomicRegion region {"X", 2000000, 2010000};
    
    //VcfReader reader {sample_vcf};
    VcfReader reader {"/Users/danielcooke/Genomics/sample_vcf/CEU.low_coverage.2010_07.xchr.genotypes.vcf"};
    
    auto variants = fetch_variants(region, reader);
    
    for (const auto& variant : variants) {
        cout << variant << endl;
    }
}

void biggest_repeats_in_huamn()
{
    auto reference = make_reference(human_reference_fasta);
    
    auto contigs = reference.get_contig_names();
    
    std::vector<GenomicRegion> biggest_regions {};
    
    auto start = std::chrono::system_clock::now();
    
    for (const auto& contig : contigs) {
        cout << "finding repeats in chromosome " << contig << endl;
        
        auto region = reference.get_contig_region(contig);
        auto repeats = find_exact_tandem_repeats(reference.get_sequence(region), region, 100);
        
        if (!repeats.empty()) {
            auto it = std::max_element(repeats.cbegin(), repeats.cend(), [] (const auto& lhs, const auto& rhs) { return size(lhs.region) < size(rhs.region); });
            biggest_regions.push_back(it->region);
        }
    }
    
    std::sort(biggest_regions.begin(), biggest_regions.end(), [] (const auto& lhs, const auto& rhs) { return size(lhs) > size(rhs); });
    
    cout << "biggest exact repeats in each chromosome" << endl;
    for (auto& region : biggest_regions) {
        cout << region << " " << size(region) << endl;
    }
    
    auto end = std::chrono::system_clock::now();
    
    auto duration = std::chrono::duration_cast<std::chrono::minutes>(end - start).count();
    
    cout << "took " << duration << " minutes" << endl;
}

void test6()
{
    auto reference = make_reference(human_reference_fasta);
    
    auto region = parse_region("Y", reference);
    
    cout << "finding repeats in " << size(region) << "bp" << endl;
    
    auto start = std::chrono::system_clock::now();
    
    auto repeats = find_exact_tandem_repeats(reference.get_sequence(region), region, 500);
    
    auto end = std::chrono::system_clock::now();
    
    auto duration = std::chrono::duration_cast<std::chrono::seconds>(end - start).count();
    
    cout << "found " << repeats.size() << " repeats" << endl;
    cout << "took " << duration << " seconds" << endl;
    
    //cout << repeats.front().region << endl;
    //cout << reference.get_sequence(repeats.front().region) << " " << repeats.front().period << endl;
    
//    auto it = std::adjacent_find(repeats.cbegin(), repeats.cend(), [] (const auto& lhs, const auto& rhs) { return lhs.region.get_end() + 5 >= rhs.region.get_begin(); });
//    
//    if (it != repeats.cend()) {
//        cout << it->region << " " << reference.get_sequence(it->region) << " " << it->period << endl;
//        ++it;
//        cout << it->region << " " << reference.get_sequence(it->region) << " " << it->period << endl;
//        ++it;
//        cout << it->region << " " << reference.get_sequence(it->region) << " " << it->period << endl;
//    }
    
//    auto it = std::adjacent_find(repeats.cbegin(), repeats.cend(), [] (const auto& lhs, const auto& rhs) { return lhs.region.get_begin() == rhs.region.get_begin(); });
//    
//    if (it != repeats.cend()) {
//        cout << it->region << " " << reference.get_sequence(it->region) << " " << it->period << endl;
//        ++it;
//        cout << it->region << " " << reference.get_sequence(it->region) << " " << it->period << endl;
//        ++it;
//        cout << it->region << " " << reference.get_sequence(it->region) << " " << it->period << endl;
//    }
    
    if (!repeats.empty()) {
        auto it = std::max_element(repeats.begin(), repeats.end(), [] (const auto& lhs, const auto& rhs) {
            return size(lhs.region) < size(rhs.region);
        });
        cout << it->region << " " << reference.get_sequence(it->region) << " " << it->period << endl;
    }
    
//    for (auto r : repeats) {
//        cout << reference.get_sequence(r) << endl;
//    }
}

int main(int argc, const char **argv)
{
    test6();
    
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
