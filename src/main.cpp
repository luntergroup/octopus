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

#include "program_options.h"
#include "octopus.h"

#include <iostream>
#include <string>
#include <chrono>
#include <algorithm>

#include "test_common.h"
#include "genomic_region.h"
#include "variant.h"
#include "reference_genome.h"
#include "read_manager.h"
#include "mappable_set.h"
#include "mappable_algorithms.h"
#include "candidate_generators.h"
#include "mappable_set.h"
#include "mappable_map.h"
#include "haplotype_tree.h"
#include "genotype_model.h"
#include "population_genotype_model.h"
#include "vcf.h"
#include "maths.h"
#include "sequence_utils.h"

#include "mock_options.h"

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
