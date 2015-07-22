//
//  input_tests.cpp
//  Octopus
//
//  Created by Daniel Cooke on 29/04/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "catch.hpp"

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <sstream>
#include <boost/program_options.hpp>

#include "test_common.h"
#include "genomic_region.h"
#include "reference_genome_factory.h"
#include "reference_genome.h"
#include "read_manager.h"
#include "program_options.h"

namespace po = boost::program_options;

using std::cout;
using std::endl;

TEST_CASE("get_reference returns a ReferenceGenome", "[program_options]")
{
    const char *argv[] = {"octopus", "--reference", human_reference_fasta.c_str(), "--reads", "test", nullptr};
    int argc = sizeof(argv) / sizeof(char*) - 1;
    
    auto options = Octopus::parse_options(argc, argv).first;
    
    auto reference = Octopus::get_reference(options);
    
    REQUIRE(reference.get_name() == "human_g1k_v37");
}

TEST_CASE("get_search_regions returns all chromosome regions when no region option is given", "[program_options]")
{
    const char *argv[] = {"octopus", "--reference", human_reference_fasta.c_str(), "--reads", "test", nullptr};
    int argc = sizeof(argv) / sizeof(char*) - 1;
    
    auto options = Octopus::parse_options(argc, argv).first;
    
    auto reference = Octopus::get_reference(options);
    
    auto regions = Octopus::get_search_regions(options, reference);
    
    auto all_contig_names = reference.get_contig_names();
    
    bool all_good = std::all_of(std::cbegin(all_contig_names), std::cend(all_contig_names),
                                [&reference, &regions] (const auto& contig_name) {
                                    return regions.count(contig_name) == 1 && reference.get_contig_region(contig_name) == *regions.at(contig_name).cbegin();
                                });
    
    REQUIRE(all_good);
}

TEST_CASE("parse_search_region_option parses manually entered regions", "[program_options]")
{
    const char *argv[] = {
        "octopus",
        "--reference", human_reference_fasta.c_str(),
        "--reads", "test",
        "--regions", "1:32000092-33000000", "5:1104209-2104209", "X:80000-900000",
        nullptr};
    int argc = sizeof(argv) / sizeof(char*) - 1;
    
    auto options = Octopus::parse_options(argc, argv).first;
    
    auto reference = Octopus::get_reference(options);
    
    auto regions = Octopus::get_search_regions(options, reference);
    
    REQUIRE(*regions["1"].cbegin() == parse_region("1:32000092-33000000", reference));
    REQUIRE(*regions["5"].cbegin() == parse_region("5:1104209-2104209", reference));
    REQUIRE(*regions["X"].cbegin() == parse_region("X:80000-900000", reference));
}

TEST_CASE("parse_search_region_option extracts regions from text files", "[program_options]")
{
    const char *argv[] = {
        "octopus",
        "--reference", human_reference_fasta.c_str(),
        "--reads", "test",
        "--regions-file", regions_txt_file.c_str(),
        nullptr};
    int argc = sizeof(argv) / sizeof(char*) - 1;
    
    auto options = Octopus::parse_options(argc, argv).first;
    
    auto reference = Octopus::get_reference(options);
    
    auto regions = Octopus::get_search_regions(options, reference);
    
    REQUIRE(regions.size() == 8); // 8 contigs
}

TEST_CASE("parse_search_region_option extracts regions from bed files", "[program_options]")
{
    const char *argv[] = {
        "octopus",
        "--reference", human_reference_fasta.c_str(),
        "--reads", "test",
        "--regions-file", regions_bed_file.c_str(),
        nullptr};
    int argc = sizeof(argv) / sizeof(char*) - 1;
    
    auto options = Octopus::parse_options(argc, argv).first;
    
    auto reference = Octopus::get_reference(options);
    
    auto regions = Octopus::get_search_regions(options, reference);
    
    REQUIRE(regions.size() == 8); // 8 contigs
}

//TEST_CASE("get_read_paths returns all read file paths", "[program_options]")
//{
//    const char *argv[] = {
//        "octopus",
//        "--reference", human_reference_fasta.c_str(),
//        "--reads", human_1000g_bam1.c_str(), human_1000g_bam2.c_str(), human_1000g_cram.c_str(),
//        nullptr};
//    int argc = sizeof(argv) / sizeof(char*) - 1;
//    
//    auto options = Octopus::parse_options(argc, argv).first;
//    
//    auto read_paths = Octopus::get_read_paths(options);
//    
//    for (auto& path : read_paths) cout << path << endl;
//}

TEST_CASE("get_read_manager returns a ReadManager with all input files", "[program_options]")
{
    const char *argv[] = {
        "octopus",
        "--reference", human_reference_fasta.c_str(),
        "--reads", human_1000g_bam1.c_str(), human_1000g_bam2.c_str(), human_1000g_cram.c_str(),
        "--reads-file", reads_file.c_str(),
        nullptr};
    int argc = sizeof(argv) / sizeof(char*) - 1;
    
    auto options = Octopus::parse_options(argc, argv).first;
    
    auto read_manager = Octopus::get_read_manager(options);
    
    REQUIRE(read_manager.get_num_samples() == 3); // the cram is same as one of the bams
}
