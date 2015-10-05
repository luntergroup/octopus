//
//  mock_options.h
//  Octopus
//
//  Created by Daniel Cooke on 14/09/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#ifndef Octopus_mock_options_h
#define Octopus_mock_options_h

#include <boost/program_options.hpp>

#include "program_options.hpp"
#include "test_common.hpp"

namespace po = boost::program_options;

inline po::variables_map get_basic_mock_options()
{
    const char *argv[] = {"octopus",
        "--reference", human_reference_fasta.c_str(),
        "--reads", human_1000g_bam1.c_str(), human_1000g_bam2.c_str(), human_1000g_bam3.c_str(),
        "--model", "population",
        //"--ploidy", "2",
        "--make-blocked-refcalls",
        //"--make-positional-refcalls",
        //"--regions", "5:157,031,410-157,031,449",
        //"--regions", "11:67503118-67503253",
        //"--regions", "2:104142870-104142984",
        //"--regions", "2:104,142,854-104,142,925",
        //"--regions", "2:104,142,897-104,142,936",
        //"--regions", "2:142376817-142376922",
        //"--regions", "20",
        //"--regions", "8:94,786,581-94,786,835", // empty region
        //"--regions", "5:157,031,410-157,031,449", "8:94,786,581-94,786,835",
        //"--regions", "5:157,031,227-157,031,282",
        //"--regions", "5:157,031,211-157,031,269",
        //"--regions", "5:157,030,955-157,031,902",
        //"--regions", "5:157,031,214-157,031,280",
        //"--regions", "2:104,142,525-104,142,742",
        "--regions", "2:104,142,403-104,142,527",
        "--min-variant-posterior", "1",
        "--min-refcall-posterior", "10",
        "--output", test_out_vcf.c_str(),
        //"--min-mapping-quality", "20",
        //"--min-snp-base-quality", "15",
        //"--no-duplicates",
        nullptr};
    
    int argc = sizeof(argv) / sizeof(char*) - 1;
    
    return Octopus::Options::parse_options(argc, argv).first;
}

#endif
