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

#include <config/common.hpp>
#include <config/option_parser.hpp>
#include <config/option_collation.hpp>
#include <logging/logging.hpp>
#include <core/octopus.hpp>

using namespace octopus;
using namespace octopus::options;

int main(const int argc, const char** argv)
{
    OptionMap options;
    
    try {
        options = parse_options(argc, argv);
    } catch (const std::exception& e) {
        logging::init();
        
        log_program_startup();
        
        logging::ErrorLogger log {};
        
        log << "A user input error has occured:";
        
        log_empty_line(log);
        
        stream(log) << '\t' << e.what();
        
        log_empty_line(log);
        
        log << "Use the command --help for descriptions of required and allowable options.";
        
        log_program_end(log);
        
        return EXIT_FAILURE;
    }
    
    if (is_run_command(options)) {
        try {
            logging::init(get_debug_log_file_name(options), get_trace_log_file_name(options));
            
            log_program_startup();
            
            run_octopus(options);
            
            log_program_end();
        } catch (const std::exception& e) {
            logging::ErrorLogger log {};
            
            log << "A program error has occured:";
            
            log_empty_line(log);
            
            stream(log) << '\t' << e.what();
            
            log_program_end(log);
            
            return EXIT_FAILURE;
        } catch (...) {
            logging::FatalLogger log {};
            log << "An unknown fatal error has occured!";
            return EXIT_FAILURE;
        }
    }
    
    return EXIT_SUCCESS;
}
