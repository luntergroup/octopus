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
#include <chrono>

#include "program_options.hpp"
#include "octopus.hpp"
#include "timing.hpp"

#include "mock_options.hpp"

int main(int argc, const char **argv)
{
    using std::cout; using std::cerr; using std::endl;
    
    auto options = get_basic_mock_options();
    
    // try {
    
    if (options) {
        if (Octopus::Options::is_run_command(*options)) {
            const auto start = std::chrono::system_clock::now();
            
            cout << "Octopus: started run at " << start << endl;
            
            Octopus::run_octopus(*options);
            
            const auto end = std::chrono::system_clock::now();
            
            cout << "Octopus: finished run at " << end << ". "
            << "Took " << TimeInterval {start, end} << endl;
        }
    } else {
        cout << "Octopus: could not parse input options. Did not start run." << endl;
    }
//    } catch (std::exception& e) {
//        cerr << "Error: " << e.what() << endl;
//        return EXIT_FAILURE;
//    } catch (...) {
//        cerr << "Error: encountered unknown error. Quiting now" << endl;
//        return EXIT_FAILURE;
//    }
    
    return EXIT_SUCCESS;
}

#endif
