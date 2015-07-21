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
#include "program_options.h"

namespace po = boost::program_options;

using std::cout;
using std::endl;

void create_argv(const std::string& args, std::vector<char*>& argv)
{
    std::istringstream iss(args);
    std::string token;
    while (iss >> token) {
        char *arg = new char[token.size() + 1];
        std::copy(token.begin(), token.end(), arg);
        arg[token.size()] = '\0';
        argv.push_back(arg);
    }
    argv.push_back(0);
}

void destroy_argv(std::vector<char*>& argv)
{
    for(std::size_t i {}; i < argv.size(); ++i) delete[] argv[i];
}

//TEST_CASE("get_search_regions returns all chromosome regions when no region option is given", "[program_options]")
//{
//    ReferenceGenomeFactory a_factory {};
//    ReferenceGenome human {a_factory.make(human_reference_fasta)};
//    
//    po::variables_map vm;
//    
//    auto regions = Octopus::get_search_regions(vm, human);
//    
//    auto all_contig_names = human.get_contig_names();
//    
//    bool all_good = std::all_of(std::cbegin(all_contig_names), std::cend(all_contig_names),
//                                [&human, &regions] (const auto& contig_name) {
//                                    return regions.count(contig_name) == 1 && human.get_contig_region(contig_name) == regions.at(contig_name).front();
//                                });
//    
//    REQUIRE(all_good);
//}

TEST_CASE("parse_search_region_option parses manually entered regions", "[program_options]")
{
    ReferenceGenomeFactory a_factory {};
    ReferenceGenome human {a_factory.make(human_reference_fasta)};
    
    std::string region_options {"-R 1:32000092-33000000 5:1104209-2104209 X:80000-900000"};
    
    std::vector<char*> argv {};
    create_argv(region_options, argv);
    
//    auto vmp = parse_options(argv.size(), &argv[0]);
//    
//    REQUIRE(vmp.second);
//    
//    auto regions = parse_region_option(vmp.first, "regions", human);
//    
//    destroy_argv(argv);
}

TEST_CASE("parse_search_region_option extracts regions from text files", "[program_options]")
{
    ReferenceGenomeFactory a_factory {};
    ReferenceGenome human {a_factory.make(human_reference_fasta)};
    
    //    po::variables_map vm;
    //    vm.insert({"regions", po::variable_value(input_regions)});
    //    po::notify(vm);
    //
    //    auto regions = parse_region_option(vm, "regions", human);
}

TEST_CASE("parse_search_region_option extracts regions from bed files", "[program_options]")
{
    ReferenceGenomeFactory a_factory {};
    ReferenceGenome human {a_factory.make(human_reference_fasta)};
    
    //    po::variables_map vm;
    //    vm.insert({"regions", po::variable_value(input_regions)});
    //    po::notify(vm);
    //
    //    auto regions = parse_region_option(vm, "regions", human);
}

TEST_CASE("parse_search_region_option works with mixture of files and manual regions", "[program_options]")
{
    ReferenceGenomeFactory a_factory {};
    ReferenceGenome human {a_factory.make(human_reference_fasta)};
    
    //    po::variables_map vm;
    //    vm.insert({"regions", po::variable_value(input_regions)});
    //    po::notify(vm);
    //
    //    auto regions = parse_region_option(vm, "regions", human);
}
