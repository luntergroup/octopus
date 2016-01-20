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
#include <chrono>
#include <utility>

#include "program_options.hpp"
#include "octopus.hpp"
#include "timing.hpp"

#include "mock_options.hpp"

int main(int argc, const char **argv)
{
    using std::cout; using std::cerr; using std::endl;
    
    try {
        //const auto options = Octopus::parse_options(argc, argv);
        const auto options = get_basic_mock_options();
        
        if (options) {
            auto start = std::chrono::system_clock::now();
            
            cout << "Octopus: started run at " << start << endl;
            
            Octopus::run_octopus(*options);
            
            auto end = std::chrono::system_clock::now();
            
            cout << "Octopus: finished run at " << end << ". ";
            cout << "Took " << TimeInterval {start, end} << endl;
        } else {
            cout << "Octopus: could not parse input options. Did not start run." << endl;
        }
    } catch (std::runtime_error& e) {
        cerr << "Error: " << e.what() << endl;
        return EXIT_FAILURE;
    }
//    catch (...) {
//        cerr << "Error: encountered unknown error. Quiting now" << endl;
//        return EXIT_FAILURE;
//    }
    
    return EXIT_SUCCESS;
}
