//
//  program_options.h
//  Octopus
//
//  Created by Daniel Cooke on 27/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__program_options__
#define __Octopus__program_options__

#include <string>
#include <vector>
#include <boost/program_options.hpp>

class ReferenceGenome;
class GenomicRegion;

namespace po = boost::program_options;

std::pair<po::variables_map, bool> parse_options(int argc, const char** argv);

std::vector<GenomicRegion> parse_region_option(const po::variables_map& options, const std::string& region_option,
                                               const ReferenceGenome& the_reference);

#endif /* defined(__Octopus__program_options__) */
