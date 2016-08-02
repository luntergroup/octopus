//
//  test_main.cpp
//  Octopus
//
//  Created by Daniel Cooke on 29/07/2016.
//  Copyright © 2016 Oxford University. All rights reserved.
//

#include <iostream>
#include <cstdlib>
#include <stdexcept>

#include "option_parser.hpp"
#include "option_collation.hpp"
#include "octopus.hpp"
#include <logging/logging.hpp>
#include "mock_options.hpp"

using namespace octopus;

int main(const int argc, const char **argv)
{
    auto options = test::get_basic_mock_options();
    
    if (options) {
        if (options::is_run_command(*options)) {
            logging::init(options::get_debug_log_file_name(*options),
                          options::get_trace_log_file_name(*options));
            run_octopus(*options);
        }
    } else {
        std::clog << "Could not parse input options. Did not run Octopus." << std::endl;
    }
}
