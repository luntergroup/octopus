//
//  program_options.h
//  Octopus
//
//  Created by Daniel Cooke on 27/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__program_options__
#define __Octopus__program_options__

#include <iostream>
#include <boost/program_options.hpp>

namespace po = boost::program_options;

std::pair<po::variables_map, bool> parse_options(int argc, char** argv);

#endif /* defined(__Octopus__program_options__) */
