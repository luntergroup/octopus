//
//  misc_benchmarks.cpp
//  Octopus
//
//  Created by Daniel Cooke on 30/08/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#define BOOST_TEST_DYN_LINK

#include <boost/test/unit_test.hpp>

#include <iostream>
#include <string>

#include "test_common.h"
#include "benchmark_utils.h"
#include "reference_genome.h"
#include "read_manager.h"
#include "variant.h"
#include "variant_utils.h"
#include "candidate_variant_generator.h"
#include "alignment_candidate_variant_generator.h"
#include "haplotype.h"
#include "genotype.h"
#include "allele.h"
#include "haplotype_tree.h"

//BOOST_AUTO_TEST_CASE(computing_hardy_weinberg)
//{
//    
//}
