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
#include "genomic_region.h"
#include "reference_genome.h"
#include "reference_genome_factory.h"
#include "read_manager.h"

using std::cout;
using std::endl;

namespace Octopus
{
    void run_octopus(po::variables_map& options)
    {
        cout << "starting Octopus" << endl;
        
//        const auto& reference_path = options.at("reference").as<std::string>();
//        
//        ReferenceGenomeFactory a_factory {};
//        
//        ReferenceGenome reference {a_factory.make(reference_path)};
//        
//        ReadManager read_manager {Octopus::get_read_paths(options)};
//        
//        SearchRegions search_regions {Octopus::get_search_regions(options, reference)};
//        
//        auto samples = read_manager.get_sample_ids();
        
    }
} // end namespace Octopus


