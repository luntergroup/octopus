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

#include <iostream>
#include <cstdlib>
#include <stdexcept>

#include "program_options.hpp"
#include "octopus.hpp"

#include <iostream>
#include <string>
#include <chrono>
#include <algorithm>

#include "test_common.hpp"
#include "genomic_region.hpp"
#include "variant.hpp"
#include "reference_genome.hpp"
#include "read_manager.hpp"
#include "mappable_set.hpp"
#include "mappable_algorithms.hpp"
#include "candidate_generators.hpp"
#include "mappable_set.hpp"
#include "mappable_map.hpp"
#include "haplotype_tree.hpp"
#include "genotype_model.hpp"
#include "population_genotype_model.hpp"
#include "vcf.hpp"
#include "maths.hpp"
#include "sequence_utils.hpp"

#include "mock_options.hpp"

using std::cout;
using std::endl;

int main(int argc, const char **argv)
{
    try {
        //auto options = Octopus::parse_options(argc, argv);
        auto options = std::make_pair(get_basic_mock_options(), true);
        
        if (options.second) {
            Octopus::run_octopus(options.first);
            std::cout << "finished running Octopus" << std::endl;
        } else {
            std::cout << "did not run Octopus" << std::endl;
        }
        
    } catch (std::runtime_error& e) {
        std::cerr << "Error: " << e.what() << std::endl;
        return EXIT_FAILURE;
    } catch (...) {
        std::cerr << "Error: encountered unknown error. Quiting now" << std::endl;
        return EXIT_FAILURE;
    }
    
    return EXIT_SUCCESS;
}
