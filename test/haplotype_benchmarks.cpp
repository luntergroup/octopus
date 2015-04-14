//
//  haplotype_benchmarks.cpp
//  Octopus
//
//  Created by Daniel Cooke on 13/04/2015.
//  Copyright (c) 2015 Oxford University. All rights reserved.
//

#include "catch.hpp"

#include <iostream>
#include <string>

#include "test_common.h"
#include "benchmark_utils.h"
#include "reference_genome.h"
#include "reference_genome_factory.h"
#include "read_manager.h"
#include "variant.h"
#include "variant_utils.h"
#include "variant_candidate_generator.h"
#include "alignment_candidate_variant_generator.h"
#include "haplotype.h"
#include "genotype.h"

//TEST_CASE("haplotype_get_sequence_benchmark", "[haplotype]")
//{
//    ReferenceGenomeFactory a_factory {};
//    ReferenceGenome human(a_factory.make(human_reference_fasta));
//    
//    auto a_region = parse_region("16:9300000-9300100", human);
//    
//    Allele allele1 {parse_region("16:9300037-9300037", human), "TG"};
//    Allele allele2 {parse_region("16:9300037-9300051", human), ""};
//    Allele allele3 {parse_region("16:9300039-9300051", human), ""};
//    
//    Haplotype hap1 {human, a_region};
//    hap1.push_back(allele3);
//    
//    Haplotype hap2 {human, a_region};
//    hap2.push_back(allele1);
//    hap2.push_back(allele2);
//    
//    auto f_get_sequence = [&]
//}
