//
//  octopus.cpp
//  Octopus
//
//  Created by Daniel Cooke on 08/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "octopus.h"

#include <iostream>

#include "common.h"
#include "program_options.h"
#include "reference_genome.h"
#include "read_manager.h"

using std::cout;
using std::endl;

namespace Octopus
{
    void run_octopus(po::variables_map& options)
    {
        cout << "starting Octopus" << endl;
        
        auto num_threads  = Octopus::get_num_threads(options);
        
        auto reference    = Octopus::get_reference(options);
        auto read_manager = Octopus::get_read_manager(options);
        auto regions      = Octopus::get_search_regions(options, reference);
        
        auto samples = read_manager.get_sample_ids();
        
    }
    
    namespace detail
    {
        void process_contig(ReferenceGenome& reference)
        {
            
        }
    }
} // end namespace Octopus


