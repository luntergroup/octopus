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
#include <boost/program_options.hpp>

#include "test_common.h"
#include "genomic_region.h"
#include "reference_genome_factory.h"
#include "reference_genome.h"
#include "program_options.h"

namespace po = boost::program_options;

using std::cout;
using std::endl;

TEST_CASE("parse_search_region_option returns all chromosome regions when no region option is given", "[program_options]")
{
    ReferenceGenomeFactory a_factory {};
    ReferenceGenome human {a_factory.make(human_reference_fasta)};
    
    po::variables_map vm;
    
    auto regions = parse_region_option(vm, "regions", human);
    
    bool has_all_chromosomes {true};
    
    for (const auto& chromosome : human.get_contig_names()) {
        auto region = human.get_contig_region(chromosome);
        
        if (std::find(regions.cbegin(), regions.cend(), region) == regions.cend()) {
            has_all_chromosomes = false;
        }
    }
    
    REQUIRE(has_all_chromosomes);
}

TEST_CASE("parse_search_region_option parses manually entered regions", "[program_options]")
{
    ReferenceGenomeFactory a_factory {};
    ReferenceGenome human {a_factory.make(human_reference_fasta)};
    
    std::string region_options {"1:32000092-33000000 5:1104209-2104209 X:80000-900000"};
    
//    const char** options = new char[][];
//    
//    options[0] = region_options.c_str();
//    
//    auto vmp = parse_options(1, options);
//    
//    REQUIRE(vmp.second);
//    
//    auto regions = parse_region_option(vmp.first, "regions", human);
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
