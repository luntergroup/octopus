//
//  octopus.cpp
//  Octopus
//
//  Created by Daniel Cooke on 08/03/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "octopus.h"

#include <iostream>
#include <thread>
#include <future>

#include "common.h"
#include "program_options.h"
#include "reference_genome.h"
#include "read_manager.h"
#include "candidate_generators.h"
#include "vcf.h"

using std::cout;
using std::endl;

namespace Octopus
{
    void run_octopus(po::variables_map& options)
    {
        cout << "starting Octopus" << endl;
        
        //auto num_system_threads = std::thread::hardware_concurrency(); // just a hint
        //if (num_system_threads == 0) num_system_threads = 1;
        
        //auto max_threads  = Octopus::get_num_threads(options);
        
        auto reference    = Octopus::get_reference(options);
        auto read_manager = Octopus::get_read_manager(options);
        auto regions      = Octopus::get_search_regions(options, reference);
        auto vcf          = Octopus::get_output_vcf(options);
        
        auto samples = read_manager.get_sample_ids();
        
    }
    
} // end namespace Octopus


