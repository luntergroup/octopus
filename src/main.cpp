//
//  main.cpp
//  Octopus
//
//  Created by Daniel Cooke on 03/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

//#define RUN_UNIT_TESTS

#ifdef RUN_UNIT_TESTS

#define BOOST_TEST_DYN_LINK
#define BOOST_TEST_MODULE Main
#include <boost/test/unit_test.hpp>

#else

#include <iostream>
#include <cstdlib>
#include <stdexcept>

#include "program_options.hpp"
#include "octopus.hpp"

#include "logging.hpp"

#include "mock_options.hpp"

int main(int argc, const char **argv)
{
    auto options = get_basic_mock_options();
    
    try {
        if (options) {
            if (Octopus::Options::is_run_command(*options)) {
                Octopus::Logging::init(Octopus::Options::get_debug_log_file_name(*options),
                                       Octopus::Options::get_trace_log_file_name(*options));
                Octopus::run_octopus(*options);
            }
        } else {
            std::clog << "Could not parse input options. Did not run Octopus." << std::endl;
        }
    } catch (const std::exception& e) {
        std::cerr << "Error: '" << e.what() << "'" << std::endl;
        return EXIT_FAILURE;
    } catch (...) {
        std::cerr << "Error: encountered unknown error. Quiting now" << std::endl;
        return EXIT_FAILURE;
    }
    
    return EXIT_SUCCESS;
}

#endif
