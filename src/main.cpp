//
//  main.cpp
//  Octopus
//
//  Created by Daniel Cooke on 03/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

//#define CATCH_CONFIG_MAIN
#define CATCH_CONFIG_RUNNER
#include "catch.hpp"

#include <iostream>

#include "program_options.h"
#include "octopus.h"

int main(int argc, char **argv)
{
    auto options = Octopus::parse_options(argc, argv);
    
    if (options.second) {
        Octopus::run_octopus(options.first);
        std::cout << "finished running Octopus" << std::endl;
    } else {
        std::cout << "did not run Octopus" << std::endl;
    }
    
    return 0;
}
