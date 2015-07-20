//
//  main.cpp
//  Octopus
//
//  Created by Daniel Cooke on 03/02/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#define CATCH_CONFIG_MAIN
//#define CATCH_CONFIG_RUNNER
#include "catch.hpp"

#include <iostream>
#include <string>
#include <iterator>
#include <vector>

#include "test_common.h"
#include "genomic_region.h"
#include "htslib_read_facade.h"
#include "read_manager.h"
#include "mock_objects.h"
#include "mappable_algorithms.h"

using std::cout;
using std::endl;

//int main(int argc, char **argv)
//{
//    ReadManager a_read_manager {std::vector<std::string> {donna_bam1, donna_bam2, donna_bam3, donna_bam4}};
//    
//    cout << a_read_manager.get_num_samples() << std::endl;
//    
//    for (const auto& sample : a_read_manager.get_sample_ids()) {
//        cout << sample << endl;
//    }
//    
//    GenomicRegion a_region {"flattened_line_28", 2000, 3000};
//    
//    cout << a_read_manager.fetch_reads(a_read_manager.get_sample_ids(), a_region).size() << endl;
//    
//    return 0;
//}
