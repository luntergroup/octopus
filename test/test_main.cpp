//
//  test_main.cpp
//  octopus
//
//  Created by Daniel Cooke on 29/07/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#include <iostream>
#include <cstdlib>
#include <stdexcept>

#include <config/common.hpp>
#include <config/option_parser.hpp>
#include <config/option_collation.hpp>
#include <core/octopus.hpp>
#include <logging/logging.hpp>

#include "mock_options.hpp"

using namespace octopus;
using namespace octopus::options;

int main(const int argc, const char** argv)
{
    OptionMap options;
    
    try {
        options = test::get_basic_mock_options();
    } catch (const std::exception& e) {
        logging::init();
        
        log_program_startup();
        
        logging::ErrorLogger error_log {};
        
        error_log << "A user input error has occured:";
        
        log_empty_line(error_log);
        
        stream(error_log) << '\t' << e.what();
        
        log_empty_line(error_log);
        
        error_log << "Use the command --help for descriptions of required and allowable options.";
        
        log_program_end(error_log);
        
        return EXIT_FAILURE;
    }
    
    if (is_run_command(options)) {
        logging::init(get_debug_log_file_name(options), get_trace_log_file_name(options));
        
        log_program_startup();
        
        run_octopus(options);
        
        log_program_end();
    }
    
    return EXIT_SUCCESS;
}
