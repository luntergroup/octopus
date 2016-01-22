//
//  octopus.hpp
//  Octopus
//
//  Created by Daniel Cooke on 08/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__octopus__
#define __Octopus__octopus__

#include <boost/program_options.hpp>

namespace po = boost::program_options;

namespace Octopus
{
    void run_octopus(po::variables_map& options);
}

#endif /* defined(__Octopus__octopus__) */
