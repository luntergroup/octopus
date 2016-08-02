//
//  main.cpp
//  Octopus
//
//  Created by Daniel Cooke on 03/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include <iostream>
#include <cstdlib>
#include <stdexcept>

#include "option_parser.hpp"
#include "option_collation.hpp"
#include "octopus.hpp"
#include <logging/logging.hpp>

using namespace octopus;

//int main(const int argc, const char **argv)
//{
//    auto options = options::parse_options(argc, argv);
//    
//    try {
//        if (options) {
//            if (options::is_run_command(*options)) {
//                logging::init(options::get_debug_log_file_name(*options),
//                              options::get_trace_log_file_name(*options));
//                run_octopus(*options);
//            }
//        } else {
//            std::clog << "Could not parse input options. Did not run Octopus." << std::endl;
//        }
//    } catch (const std::exception& e) {
//        std::cerr << "Error: '" << e.what() << "'" << std::endl;
//        return EXIT_FAILURE;
//    } catch (...) {
//        std::cerr << "Encountered unknown error. Quiting now" << std::endl;
//        return EXIT_FAILURE;
//    }
//    
//    return EXIT_SUCCESS;
//}
