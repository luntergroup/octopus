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
#include <memory>

#include "common.h"
#include "program_options.h"
#include "mappable_map.h"
#include "reference_genome.h"
#include "read_manager.h"
#include "read_filter.h"
#include "read_transform.h"
#include "read_utils.h"
#include "candidate_generators.h"
#include "vcf.h"
#include "variant_caller.h"
#include "basic_caller.h"

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
        
        auto reference           = Octopus::get_reference(options);
        auto read_manager        = Octopus::get_read_manager(options);
        auto regions             = Octopus::get_search_regions(options, reference);
        auto read_filter         = Octopus::get_read_filter(options);
        auto read_transform      = Octopus::get_read_transformer(options);
        auto candidate_generator = Octopus::get_candidate_generator(options, reference);
        auto vcf                 = Octopus::get_output_vcf(options);
        
        for (const auto& contig_region : regions) {
            auto region = *contig_region.second.cbegin();
            
            cout << "looking at region " << region << endl;
            
            std::unique_ptr<VariantCaller> caller = std::make_unique<BasicVariantCaller>(reference, read_manager, read_filter, read_transform, candidate_generator);
            
            auto calls = caller->call_variants(region);
        }
    }
    
} // end namespace Octopus


