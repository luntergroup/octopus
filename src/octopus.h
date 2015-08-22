//
//  octopus.h
//  Octopus
//
//  Created by Daniel Cooke on 08/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef __Octopus__octopus__
#define __Octopus__octopus__

#include <vector>
#include <boost/program_options.hpp>

#include "genomic_region.h"

namespace po = boost::program_options;

class ReferenceGenome;

namespace Octopus
{
    void run_octopus(po::variables_map& options);
    
    namespace detail
    {
        void process_contig(ReferenceGenome& reference);
    }
} // end namespace Octopus

#endif /* defined(__Octopus__octopus__) */
